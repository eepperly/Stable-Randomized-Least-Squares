%%%% Convergence of FOSSILS
%%%% Iteration count grid

clear,clc
addpath('../code')
addpath('../utils')

%% Input
% Input problem
cond_A = 1e8;
res_size = 1e-3;
s = 12;

M = 15; N = 21;
mvec = round(logspace(3,6,M));
nvec = round(logspace(log10(50),log10(1e3),N));

% Monitor convergence and operation count (roughly)    
itcount = zeros(M,N);

%% Loop through conds and residuals
for ii = 1:M
    m = mvec(ii);

    currentnvec = nvec(nvec <= m/s);
    nmax = currentnvec(end);
    Ufull = haarorth(m,nmax+1);
    
    'starting'
    m
    currentnvec
    for jj = 1:length(currentnvec)
        n = nvec(jj);

        % Make data
        [A,b,x,r] = random_ls_problem(m,n,cond_A,res_size,Ufull);

        [y,~,num_iters] = spir(A,b,[],[],[],[],true);
        itcount(ii,jj) = num_iters;

        itcount
    end
end
save('../data/grid_size_SPIR.mat',"mvec","nvec","cond_A","res_size","itcount")

%% Plot
[nvec_plt,mvec_plt] = meshgrid(nvec',mvec);

figure
datamat = itcount;
datamat(datamat == 0) = NaN;
contourf(nvec_plt,mvec_plt,datamat),
set(gca,'fontsize', 16)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('$n$', 'Interpreter','latex', 'fontsize', 20)
ylabel('$m$', 'Interpreter','latex', 'fontsize', 20)
title("Iteration count",'fontsize', 16,'Interpreter','latex')
colorbar, colormap('cool')

saveas(gcf,'../figs/grid_size.fig')
