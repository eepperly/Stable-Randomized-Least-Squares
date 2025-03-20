function [x,stats,num_iters] = metasolver(A,b,setup,iterate,varargin)
    Anum = isnumeric(A);
    if Anum
        scale = vecnorm(A);
        A = A ./ scale;
        m = size(A,1);
        n = size(A,2);
        Afun = @(x,op) mul(A,x,op);
    else
        scale = 1;
        m = size(b,1);
        n = size(A(zeros(size(b,1),0), true),1);
        Afun = A;
    end

    if length(varargin) >= 1 && ~isempty(varargin{1})
        d = varargin{1};
    else
        d = 12*n;
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        iterations = varargin{2};
    else
        iterations = 'adaptive';
    end
    
    if length(varargin) >= 3 && ~isempty(varargin{3})
        summary = varargin{3};
    else
        summary = [];
    end

    if length(varargin) >= 4 && ~isempty(varargin{4})
        verbose = varargin{4};
    else
        verbose = false;
    end

    if length(varargin) >= 5 && ~isempty(varargin{5})
        reproducible = varargin{5};
    else
        reproducible = false;
    end

    if length(varargin) >= 6 && ~isempty(varargin{6})
        lowrankprecon = varargin{6};
    else
        lowrankprecon = false;
    end
    
    if ~reproducible && exist('sparsesign','file') == 3
        S = sparsesign(d,m,8);
    else
        warning(['Using slower and slightly incorrect backup ' ...
            'implementation sparse_sign_backup.m. For the better ' ...
            'implementation, build the mex file `sparsesign.c` ' ...
            'using the command `mex sparsesign.c`.']);
        S = sparse_sign_backup(d,m,8);
    end

    betol = eps;
    if ~isstring(iterations) && ~ischar(iterations)
        if length(iterations) == 1 && (floor(iterations) ~= iterations)
            betol = iterations;
            iterations = [100,100];
            adaptive_switching = true;
            adaptive_stopping = true;
        else
            adaptive_switching = false;
            adaptive_stopping = false;
        end
    elseif strcmp(iterations, 'halfadaptive')
        adaptive_switching = true;
        adaptive_stopping = false;
        iterations = [100,30];
    elseif strcmp(iterations, 'adaptive')
        adaptive_switching = true;
        adaptive_stopping = true;
        iterations = [100,100];
    else
        error("Second argument should be a list of integers or 'adaptive'")
    end

    stats = [];

    if Anum
        SA = full(S*A);
    else
        SA = full(Afun(S',true)');
    end
    [U,svals,V] = svd(SA,'econ'); svals = diag(svals);

    % Norm and condition number estimation
    Acond = max(svals) / min(svals);
    Anorm = max(svals); 
    Afronorm = norm(svals);
    bnorm = norm(b);
    ftol = 10*(Anorm*norm(x) + 0.04*Acond*norm(r))*eps;

    if Acond > 1/eps/30
        if lowrankprecon
            ind = svals/max(svals) > eps*30;
            sreg = svals(ind);
            Vreg = V(:, ind);
            Ureg = U(:, ind);
        else
            reg = 10 * Afronorm * eps;
            sreg = sqrt(svals.^2 + reg^2);
            Vreg = V;
            Ureg = U;
        end
    else
        reg = 0;
        sreg = svals;
        Vreg = V;
        Ureg = U;
    end

    % Initial guess
    x = Vreg*((Ureg'*(S*b)) ./ sreg);
    r = b - Afun(x,false);

    betol = Afronorm * betol; % Backward error tolerance

    if ~isempty(summary); stats(end+1,:) = summary(x./scale.'); end
    num_iters = 0;
    stopnow = false;
    stopnext = false;

    for loop = 1:length(iterations)
        c = (Vreg'*(Afun(r,true) - reg^2*x))./sreg;
        dy = c;
        solverdata = setup(c,dy,[r;-reg*x],...
            @(xx) RAAR(xx,Vreg,sreg,Afun,reg),...
            @(xx) AR(xx,Vreg,sreg,Afun,reg),...
            @(xx) RA(xx,Vreg,sreg,Afun,reg),d);
        if verbose; fprintf('Beginning loop %d\n', loop); end
        for i = 1:iterations(loop)
            num_iters = num_iters + 1;
            [dy,update,solverdata] = iterate(solverdata);
            if ~isempty(summary)
                stats(end+1,:) = summary((x + V*(dy./sreg)) ...
                    ./ scale.'); %#ok<AGROW>
            end
            if verbose
                xhat = x + Vreg*(dy./sreg);
                rhat = b - Afun(xhat,false);
                be = posterior_estimate(Afun,xhat,rhat,V,svals,Afronorm,bnorm);
                fprintf('Iteration %d\t%e / %e\t%e / %e\n', i,...
                    norm(update), ftol, be, betol);
            end
            if adaptive_switching && loop == 1 && norm(update) <= ftol
                break
            elseif adaptive_stopping && loop == 2 && mod(i,5) == 0
                xhat = x + Vreg*(dy./sreg);
                rhat = b - Afun(xhat,false);
                be = posterior_estimate(Afun,xhat,rhat,V,svals,Afronorm,bnorm);
                if be <= betol || stopnext
                    stopnow = true; 
                    x = xhat;
                    break; 
                elseif be <= 4*betol
                    stopnext = true;
                end
            end
        end
        if stopnow; break; end
        x = x + Vreg*(dy./sreg);
        r = b - Afun(x,false);
        if adaptive_stopping
            be = posterior_estimate(Afun,x,r,V,svals,Afronorm,bnorm);
            if be <= betol; break; end
            if loop == 2
                warning('Exiting without backward stability!')
            end
        end
    end

    x = x ./ scale.';
end

function be = posterior_estimate(Afun,x,r,V,svals,Afronorm,bnorm)
    theta = Afronorm / bnorm;
    xnorm = norm(x);
    be = norm((V'*Afun(r,true)) ./ (svals.^2 ...
        + theta^2*norm(r)^2/(1+theta^2*xnorm^2)).^0.5) ...
        * (theta / sqrt(1+theta^2*xnorm^2));
end

function b = mul(A, x, op)
    if op
        b = A'*x;
    else
        b = A*x;
    end
end

function RAARdy = RAAR(dy,V,sreg,Afun,reg)
    z = V*(dy./sreg);
    RAARdy = (V'*(Afun(Afun(z,false),true) + reg^2 * z)) ./ sreg;
end

function ARdy = AR(dy,V,sreg,Afun,reg)
    z = V*(dy./sreg);
    ARdy = [Afun(z,false);reg*z];
end

function RAb = RA(b,V,sreg,Afun,reg)
    z = Afun(b(1:(end-size(V,1)),:), true)+reg*b((end-size(V,1)+1):end,:);
    RAb = (V'*z) ./sreg;
end

