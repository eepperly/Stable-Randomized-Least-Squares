%%%% Convergence of FOSSILS
%%%% Iteration count grid

clear,clc
addpath('../code')
addpath('../utils')

%% Input
% Input problem
N = 30;

conds = logspace(0,15,N);
ress = logspace(-15,0,N);

m = 4000; n = 50;

% Monitor convergence and operation count (roughly)    
bwerr_foss = zeros(N,N);
bwerr_ism = zeros(N,N);
itcount = zeros(N,N);

U = haarorth(m,n+1);
%% Loop through conds and residuals
for ii = 1:N
    cond_A = conds(ii);
    for jj = 1:N
        res_size = ress(jj);

        % Make data
        [A,b,x,r] = random_ls_problem(m,n,cond_A,res_size,U);
        
        [y,~,num_iters] = fossils(A,b);

        bwerr_foss(ii,jj) = kw_estimate(A,b,y)/norm(A,'fro');
        itcount(ii,jj) = num_iters;

        y = iterative_sketching(A,b,[],[],[],[],"optimal","optimal");
        bwerr_ism(ii,jj) = kw_estimate(A,b,y)/norm(A,'fro');

        itcount
    end
end
save('../data/grid_cond_res_FOSSILS.mat',"m","n", "ress","conds",...
    "itcount","bwerr_foss","bwerr_ism")

%% Plot

% First plot: backward error

[ress_plt,conds_plt] = meshgrid(ress',conds);
figure('Position', [353   452   774   303])

be_low = 1e-18; be_high = 1e-10;
clevels = logspace(-18,-10,15);

subplot(121)
datamat = bwerr_foss;
contourf(ress_plt,conds_plt,datamat,clevels)
colormap('cool')
set(gca,'fontsize', 14)
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
set(gca,'ColorScale','log')
colorbar('eastoutside')
clim([be_low, be_high])
xlabel('Residual $\|\mbox{\boldmath $b$} - \mbox{\boldmath $A$}\mbox{\boldmath $x$}\|$', 'Interpreter','latex', 'fontsize', 14)
ylabel('Condition number $\kappa$', 'Interpreter','latex', 'fontsize', 14)
title("FOSSILS",'fontsize', 14,'Interpreter','latex')

subplot(122)
datamat = bwerr_ism;
contourf(ress_plt,conds_plt,datamat,clevels)
colormap('cool')
set(gca,'fontsize', 14)
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
set(gca,'ColorScale','log')
colorbar('eastoutside')
clim([be_low, be_high])
xlabel('Residual $\|\mbox{\boldmath $b$} - \mbox{\boldmath $A$}\mbox{\boldmath $x$}\|$', 'Interpreter','latex', 'fontsize', 14)
ylabel('Condition number $\kappa$', 'Interpreter','latex', 'fontsize', 14)
title("Iter sketch + mom",'fontsize', 14,'Interpreter','latex')

saveas(gcf,'../figs/backerrs.fig')
saveas(gcf,'../figs/backerrs.png')

% Second plot: iterations

[ress_plt,conds_plt] = meshgrid(ress',conds);
figure('Position', [353   452   774   303])
iter_low = 5;
iter_high = 5 * ceil(max(itcount(:))/5);

subplot(121)
datamat = itcount;
contourf(ress_plt,conds_plt,datamat,15)
set(gca,'fontsize', 14)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('Residual $\|\mbox{\boldmath $b$} - \mbox{\boldmath $A$}\mbox{\boldmath $x$}\|$', 'Interpreter','latex', 'fontsize', 14)
ylabel('Condition number $\kappa$', 'Interpreter','latex', 'fontsize', 14)
title("Iterations ($m=4000,n=50$)",'fontsize', 14,'Interpreter','latex')

colormap cool
clim([iter_low iter_high])
colorbar('eastoutside')

try
    load('../data/grid_size_FOSSILS.mat')
    
    subplot(122)
    [nvec_plt,mvec_plt] = meshgrid(nvec',mvec);
    datamat = itcount;
    datamat(datamat == 0) = NaN;
    contourf(nvec_plt,mvec_plt,datamat,15)
    set(gca,'fontsize', 14)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlabel('Columns $n$', 'Interpreter','latex', 'fontsize', 14)
    ylabel('Rows $m$', 'Interpreter','latex', 'fontsize', 14)
    title("Iterations ($\kappa=10^8, \|\mbox{\boldmath $b$}-\mbox{\boldmath $Ax$}\| = 10^{-3}$)",'fontsize', 14,'Interpreter','latex')
    colormap cool
    clim([iter_low iter_high])
    colorbar('eastoutside')
catch
    warning("Run 'test_sizes' to populate the right plot of this figure")
end

saveas(gcf,'../figs/itcounts.fig')
saveas(gcf,'../figs/itcounts.png')
