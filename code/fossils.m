function [x,stats,num_iters] = fossils_op(A,b,varargin)
%FOSSILS Solve A*x = b in the least-squares sense using FOSSILS
%   Optional parameters (use [] for default value):
%   - d: sketching dimension (default 12*size(A,2))
%   - iterations: vector of integers. iterations(i) contains the number of
%     steps to be taken on the i-th refinement step. Or set to 'adaptive'
%     to automatically determine the number of refinement steps needed to
%     achieve backward stability. (default 'adaptive')
%   - summary: a function of the current iterate to be recorded at each
%     iteration. All summary values will be rows of stats, the second
%     output of this function.
%   - verbose: if true, print at each iteration. (default false)
%   - reproducible: if true, use slower, reproducible implementation of
%     sparse sign embeddings (default false)
    Anum = isnumeric(A);
    if Anum
        m = size(A,1);
        n = size(A,2);
        Afun = @(x,op) mul(A,x,op);
    else
        m = size(b,1);
        n = size(A(b, true),1);
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
    
    if ~reproducible && exist('sparsesign','file') == 3
        S = sparsesign(d,m,8);
    else
        warning(['Using slower and slightly incorrect backup ' ...
            'implementation sparse_sign_backup.m. For the better ' ...
            'implementation, build the mex file `sparsesign.c` ' ...
            'using the command `mex sparsesign.c`.']);
        S = sparse_sign_backup(d,m,8);
    end

    if ~isstring(iterations) && ~ischar(iterations)
        adaptive = false;
    elseif strcmp(iterations, 'adaptive')
        adaptive = true;
        iterations = [100,100];
    else
        error("Second argument should be a list of integers or 'adaptive'")
    end

    stats = [];

    if Anum
        SA = full(S*A);
    else
        SA = full(Afun(A',true)');
    end
    [U,svals,V] = svd(SA,'econ'); svals = diag(svals);
    x = V*((U'*(S*b)) ./ svals);
    r = b - Afun(x,false);

    % Chebyshev interpolation
    momentum = n/d;
    damping = (1 - momentum)^2;

    % Norm and condition number estimation
    Acond = max(svals) / min(svals);
    Anorm = max(svals); 
    if Anum
        Afronorm = norm(A,'fro');
    else
        Afronorm = norm(svals);
    end
    bnorm = norm(b);

    if Acond > 1/eps/30
        reg = 10 * Afronorm * eps;
        sreg = sqrt(svals.^2 + reg^2);
    else
        reg = 0;
        sreg = svals;
    end

    if ~isempty(summary); stats(end+1,:) = summary(x); end
    num_iters = 0;

    for loop = 1:length(iterations)
        c = (V'*(Afun(r,true) - reg^2*x))./sreg;
        dy = c; dyold = dy;

        if verbose; fprintf('Beginning loop %d\n', loop); end
        for i = 1:iterations(loop)
            num_iters = num_iters + 1;
            z = V*(dy./sreg);
            update = damping * (c - (V'*(Afun(Afun(z,false),true)...
                + reg^2 * z)) ./sreg) + momentum*(dy-dyold);
            dyold = dy;
            dy = dy + update;
            if ~isempty(summary)
                stats(end+1,:) = summary(x + V*(dy./sreg));
            end
            if verbose
                xhat = x + V*(dy./sreg);
                rhat = b - Afun(xhat,false);
                be = posterior_estimate(Afun,xhat,rhat,V,svals,Afronorm,bnorm);
                fprintf('Iteration %d\t%e\t%e\n', i, norm(update), be);
            end
            if adaptive && loop == 1 && norm(update) ...
                    <= 10*(Anorm*norm(x) + 0.04*Acond*norm(r))*eps
                break
            elseif adaptive && loop == 2 && mod(i,5) == 0
                xhat = x + V*(dy./sreg);
                rhat = b - A*xhat;
                be = posterior_estimate(Afun,xhat,rhat,V,svals,Afronorm,bnorm);
                if be <= Afronorm * eps/2; break; end
            end
        end
        x = x + V*(dy./sreg);
        r = b - Afun(x,false);
        if adaptive
            be = posterior_estimate(Afun,x,r,V,svals,Afronorm,bnorm);
            if be <= Afronorm * eps/2; break; end
            if loop == 2
                warning('Exiting without backward stability!')
            end
        end
    end
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