function [x,stats,num_iters] = spir(A,b,varargin)
%SPIR Solve A*x = b in the least-squares sense using SPIR
%(sketch-and-precondition with iterative refinement)
%   Optional parameters (use [] for default value):
%   - d: sketching dimension (default 12*size(A,2))
%   - iterations: can be set in three ways
%        * if iterations is a vector of integers, iterations(i) contains
%          the number of steps to be taken on the i-th refinement step.
%        * if iterations is a non-integer, iterations sets the tolerance
%          for the backward error and the number of steps is determined
%          adaptively
%        * if iterations is 'adaptive', the number of steps is determined
%          adaptively with the default tolerance (default 'adaptive')
%   - summary: a function of the current iterate to be recorded at each
%     iteration. All summary values will be rows of stats, the second
%     output of this function.
%   - solver: either 'cg' or 'lsqr' (default 'cg')
%   - verbose: if true, print at each iteration. (default false)
%   - reproducible: if true, use slower, reproducible implementation of
%     sparse sign embeddings (default false)
%   - truncprecon: if true, use a non-square factorization for the
%     preconditioner in case of severe ill-conditioning. If false, use
%     explicit regularization in the singular values. (default true)

    if length(varargin) >= 4
        if ~isempty(varargin{4})
            solver = varargin{4};
        else
            solver = 'cg';
        end
        varargin(4) = [];
    else
        solver = 'cg';
    end

    if strcmp(solver, 'cg')
        [x,stats,num_iters] = metasolver(A,b,@spir_setup,...
            @spir_iterate,varargin{:});
    elseif strcmp(solver,'lsqr')
        [x,stats,num_iters] = metasolver(A,b,@spir_lsqr_setup,...
            @spir_lsqr_iterate,varargin{:});
    else
        error("Solver '%s' not recognized", solver);
    end
end

function data = spir_setup(c,dy,r,RAAR,AR,RA,d) %#ok<INUSD>
    data.r = c - RAAR(dy);
    data.rsq = data.r' * data.r;
    data.p = data.r;
    data.dy = dy;
    data.matvec = RAAR;
end

function [dy,update,data] = spir_iterate(data)
    Ap = data.matvec(data.p);
    alpha = data.rsq / (data.p' * Ap);
    update = alpha * data.p;
    data.dy = data.dy + update;
    dy = data.dy;
    data.r = data.r - alpha * Ap;
    new_rsq = data.r'*data.r;
    beta = new_rsq / data.rsq;
    data.rsq = new_rsq;
    data.p = data.r + beta*data.p;
end

function data = spir_lsqr_setup(c,dy,r,RAAR,AR,RA,d) %#ok<INUSD>
    data.dy = dy;
    b = r - AR(dy);
    beta = norm(b); data.u = b / beta;
    data.v = RA(data.u);
    data.alpha = norm(data.v); data.v = data.v / data.alpha;
    data.w = data.v;
    data.phibar = beta; data.rhobar = data.alpha;
    data.AR = AR; data.RA = RA;
end

function [dy,update,data] = spir_lsqr_iterate(data)
    data.u = data.AR(data.v) - data.alpha*data.u;
    beta = norm(data.u); data.u = data.u / beta;
    data.v = data.RA(data.u) - beta*data.v;
    data.alpha = norm(data.v); data.v = data.v / data.alpha;
    rho = sqrt(data.rhobar^2 + beta^2);
    c = data.rhobar / rho;
    s = beta / rho;
    theta = s * data.alpha;
    data.rhobar = - c * data.alpha;
    phi = c * data.phibar;
    data.phibar = s * data.phibar;
    update = (phi/rho) * data.w;
    data.dy = data.dy + update; dy = data.dy;
    data.w = data.v - (theta/rho) * data.w;
end
