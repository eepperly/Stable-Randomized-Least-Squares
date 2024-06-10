function [A,b,x,r,S,V] = random_ls_problem(m,n,cond_A,res_size,varargin)
if ~isempty(varargin)
    U = varargin{1};
else
    U = haarorth(m,n+1);
end
V = haarorth(n);
S = diag(logspace(-log10(cond_A),0,n));
A = U(:,1:n)*S*V';
x = orth(randn(n,1));
r = U(:,end) * res_size;
b = A*x + r;
end