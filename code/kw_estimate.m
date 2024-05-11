function be = kw_estimate(A, b, y, varargin)
    if isempty(varargin)
        theta = Inf;
    else
        theta = varargin{1};
    end

    if length(varargin) >= 2 && strcmp(varargin{2}, "sketched")
        if length(varargin) >= 3
            d = varargin{3};
        else
            d = 2*size(A,2);
        end
        sketch = sparsesign(d,size(A,1),8);
        [~,S,V] = svd(full(sketch*A),'econ');
    elseif length(varargin) >= 3
        S = varargin{2};
        V = varargin{3};
    else
        [~,S,V] = svd(A,'econ');
    end
    
    r = b-A*y;
    if isinf(theta)
        omega = norm(r) / norm(y);
    else
        omega = theta*norm(r) / sqrt(1 + theta^2*norm(x)^2);
    end

    be = omega / norm(r) * norm((V'*(A'*r)) ./ (diag(S).^2 + omega^2).^0.5);
end