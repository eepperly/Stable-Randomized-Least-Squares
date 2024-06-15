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
%   - verbose: if true, print at each iteration. (default false)
%   - reproducible: if true, use slower, reproducible implementation of
%     sparse sign embeddings (default false)

    [x,stats,num_iters] = metasolver(A,b,@spir_setup,...
        @spir_iterate,varargin{:});
end

function data = spir_setup(c,dy,matvec,d) %#ok<INUSD>
    data.r = c - matvec(dy);
    data.rsq = data.r' * data.r;
    data.p = data.r;
    data.dy = dy;
    data.matvec = matvec;
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
