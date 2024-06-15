function [x,stats,num_iters] = fossils(A,b,varargin)
%FOSSILS Solve A*x = b in the least-squares sense using FOSSILS
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
    [x,stats,num_iters] = metasolver(A,b,@fossils_setup,...
        @fossils_iterate,varargin{:});
end

function data = fossils_setup(c,dy,matvec,d)
    data.dyold = dy;
    data.dy = dy;
    data.c = c;
    data.matvec = matvec;
    data.momentum = size(c,1)/d;
    data.damping = (1-data.momentum)^2;
end

function [dy,update,data] = fossils_iterate(data)
    Ady = data.matvec(data.dy);
    update = data.damping * (data.c - Ady)...
        + data.momentum * (data.dy - data.dyold);
    data.dyold = data.dy;
    data.dy = data.dy + update;
    dy = data.dy;
end

