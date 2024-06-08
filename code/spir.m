function [x,stats] = spir(A,b,varargin)
%SPIR Solve A*x = b in the least-squares sense by sketch-and-precondition
%with iterative refinement
%   Optional parameters (use [] for default value):
%   - d: sketching dimension (default 2*n).
%   - q: number of steps.
%   - summary: a function of the current iterate to be recorded at each
%     iteration. All summary values will be rows of stats, the second
%     output of this function.
%   - verbose: if true, print at each iteration.
%   - opts: specify the solver ('cgne' or 'lsqr') and whether to use a
%     'warm' start (initial iterate given by sketch-and-solve) or a 'cold'
%     start (initial iterate zero). Defaults to 'lsqrwarm'
%   - reproducible: if true, use a slow, but reproducible implementation
%     of sparse sign embeddings. Defaults to false

    m = size(A,1);
    n = size(A,2);

    if length(varargin) >= 1 && ~isempty(varargin{1})
        d = varargin{1};
    else
        d = 2*n;
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        qs = varargin{2};
    else
        qs = [50,50];
    end
    
    if length(varargin) >= 3 && ~isempty(varargin{3})
        summary = varargin{3};
    else
        summary = [];
        mysum = [];
    end

    if length(varargin) >= 4 && ~isempty(varargin{4})
        verbose = varargin{4};
    else
        verbose = false;
    end

    if length(varargin) >= 5 && ~isempty(varargin{5})
        opts = varargin{5};
    else
        opts = '';
    end

    if length(varargin) >= 6 && ~isempty(varargin{6})
        reproducible = varargin{6};
    else
        reproducible = false;
    end
    
    if ~reproducible && exist('sparsesign','file') == 3
        S = sparsesign(d,m,8);
    else
        warning(['Using slow implementation sparsesign_slow.m. ' ...
            'For the better ' ...
            'implementation, build the mex file `sparsesign.c` ' ...
            'using the command `mex sparsesign.c`.']);
        S = sparsesign_slow(d,m,8);
    end

    [Q,R] = qr(full(S*A),'econ');

    if contains(opts, 'cold')
        x = zeros(n,1);
    else
        x = R \ (Q' * (S*b));
    end

    for q_idx = 1:length(qs)
        q = qs(q_idx);
        if contains(opts, 'cgne')
            if ~isempty(summary)
                mysym = @(y) summary(x+y);
            end
            [dx,~,newstats] = mycg(@(y) A'*(A*y),@(y) R\(R'\y),...
                A'*(b-A*x),0,q,mysym,zeros(size(x)),verbose);
        elseif contains(opts, 'sym')
            if ~isempty(summary)
                mysum = @(y) summary(x+R\y);
            end
            [dy,~,newstats] = mycg(@(y) R'\(A'*(A*(R\y))),...
                @(y) y,R'\(A'*(b-A*x)),0,q,mysum,zeros(size(x)),verbose);
            dx = R\dy;
        else
            if ~isempty(summary)
                mysum = @(y) summary(x+R\y);
            end
            [dy,~,newstats] = mylsqr(@(y) A*(R\y),@(y) R'\(A'*y),b-A*x,...
                0,q,mysum,zeros(size(x)),verbose);
            dx = R\dy;
        end
        x = x + dx;
        if ~isempty(summary)
            if q_idx == 1
                stats = newstats;
            else
                stats = [stats;newstats(2:end,:)]; %#ok<AGROW>
            end
        end
    end

end