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
nvec = round([50, logspace(2,4,N-1)]);

% Plotting input
MS = 'Markersize'; LW = 'linewidth'; FS = 'fontsize';
TEX = 'interpreter';tex = 'latex'; CO = 'Color';
lw = 2; ms = 12; fs = 14; ffs = 10; lww = 3; mss = 24; fss = 20; 
co = {'b','r','k','m',[0 .6 0],[.6 0 0],[0 0 .6]}; linestyles = {'-','--',':','-.'}; markers = {'o','x','^','s'};
green = [0 .7 0];

% Monitor convergence and operation count (roughly)    
bwerr = zeros(M, N);
residuals = zeros(M, N);
itcount = zeros(M,N);
forerr = zeros(M,N);

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

        [y,~,num_iters] = fossils(A,b);

        summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);backward_error_ls(A,b,y)/norm(A,'fro')];
        sumit = summary(y);

        bwerr(ii,jj) = sumit(3);
        residuals(ii,jj) = sumit(2);
        forerr(ii,jj) = sumit(1);
        itcount(ii,jj) = num_iters;

        itcount
    end
end
save('grid_size_FOSSILS_June3_THROWAWAY.mat',"mvec","nvec", "cond_A","res_size","itcount","forerr", "residuals","bwerr")

%% Plot
[nvec_plt,mvec_plt] = meshgrid(nvec',mvec);

figure
subplot(221)
datamat = itcount;
datamat(datamat == 0) = NaN;
contourf(log10(nvec_plt),log10(mvec_plt),datamat),
set(gca,'fontsize', 16)
xlabel('$n$', 'Interpreter','latex', 'fontsize', 20)
ylabel('$m$', 'Interpreter','latex', 'fontsize', 20)
title("Iteration count",'fontsize', 16,'Interpreter','latex')
colorbar, colormap('cool')

subplot(222)
datamat = bwerr;
datamat(datamat == 0) = NaN;
contourf(log10(nvec_plt/eps),log10(mvec_plt),log10(datamat)),
set(gca,'fontsize', 16)
xlabel('$n$', 'Interpreter','latex', 'fontsize', 20)
ylabel('$m$', 'Interpreter','latex', 'fontsize', 20)
title("Backward error $\mbox{BE}_\infty(\mbox{\boldmath $\widehat{x}$}_i)$",'fontsize', 16,'Interpreter','latex')
colorbar, colormap('cool')

subplot(223)
datamat = forerr;
datamat(datamat == 0) = NaN;
contourf(log10(nvec_plt/eps),log10(mvec_plt),log10(datamat)),
set(gca,'fontsize', 16)
xlabel('$n$', 'Interpreter','latex', 'fontsize', 20)
ylabel('$m$', 'Interpreter','latex', 'fontsize', 20)
title("Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$",'fontsize', 16,'Interpreter','latex')
colorbar, colormap('cool')

subplot(224)
datamat = residuals;
datamat(datamat == 0) = NaN;
contourf(log10(nvec_plt/eps),log10(mvec_plt),log10(datamat)),
set(gca,'fontsize', 16)
xlabel('$n$', 'Interpreter','latex', 'fontsize', 20)
ylabel('$m$', 'Interpreter','latex', 'fontsize', 20)
title("Residual error $\|\mbox{\boldmath $A$}(\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i)\|/\|\mbox{\boldmath $b$}\|$",'fontsize', 16,'Interpreter','latex')
colorbar, colormap('cool')

saveas(gcf,'../figures/grid_size.fig')
