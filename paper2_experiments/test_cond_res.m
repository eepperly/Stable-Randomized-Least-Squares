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

% Plotting input
MS = 'Markersize'; LW = 'linewidth'; FS = 'fontsize';
TEX = 'interpreter';tex = 'latex'; CO = 'Color';
lw = 2; ms = 12; fs = 14; ffs = 10; lww = 3; mss = 24; fss = 20; 
co = {'b','r','k','m',[0 .6 0],[.6 0 0],[0 0 .6]}; linestyles = {'-','--',':','-.'}; markers = {'o','x','^','s'};
green = [0 .7 0];

% Monitor convergence and operation count (roughly)    
bwerr = zeros(N, N);
residuals = zeros(N, N);
itcount = zeros(N,N);
forerr = zeros(N,N);

U = haarorth(m,n+1);
%% Loop through conds and residuals
for ii = 1:N
    cond_A = conds(ii);
    for jj = 1:N
        res_size = ress(jj);

        % Make data
        [A,b,x,r] = random_ls_problem(m,n,cond_A,res_size,U);
        
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
save('grid_cond_res_FOSSILS_June3.mat',"m","n", "ress","conds","itcount","forerr", "residuals","bwerr")

%% Plot
[ress_plt,conds_plt] = meshgrid(ress',conds);

figure
subplot(221)
datamat = itcount;
contourf(log10(ress_plt),log10(conds_plt),datamat),
set(gca,'fontsize', 16)
xlabel('$\|Ax^*-b\|_2$', 'Interpreter','latex', 'fontsize', 20)
ylabel('$\kappa_2(A)$', 'Interpreter','latex', 'fontsize', 20)
title("Iteration count",'fontsize', 16,'Interpreter','latex')
colorbar, colormap('cool')

subplot(222)
datamat = bwerr;
contourf(log10(ress_plt/eps),log10(conds_plt),log10(datamat)),
set(gca,'fontsize', 16)
xlabel('$\|Ax-b\|/u$', 'Interpreter','latex', 'fontsize', 20)
ylabel('$\kappa$', 'Interpreter','latex', 'fontsize', 20)
title("Backward error $\mbox{BE}_\infty(\mbox{\boldmath $\widehat{x}$}_i)$",'fontsize', 16,'Interpreter','latex')
colorbar, colormap('cool')

subplot(223)
datamat = forerr;
datamat(datamat == 0) = NaN;
contourf(log10(ress_plt/eps),log10(conds_plt),log10(datamat)),
set(gca,'fontsize', 16)
xlabel('$\|Ax-b\|/u$', 'Interpreter','latex', 'fontsize', 20)
ylabel('$\kappa$', 'Interpreter','latex', 'fontsize', 20)
title("Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$",'fontsize', 16,'Interpreter','latex')
colorbar, colormap('cool')

subplot(224)
datamat = residuals;
contourf(log10(ress_plt/eps),log10(conds_plt),log10(datamat)),
set(gca,'fontsize', 16)
xlabel('$\|Ax-b\|/u$', 'Interpreter','latex', 'fontsize', 20)
ylabel('$\kappa$', 'Interpreter','latex', 'fontsize', 20)
title("Residual error $\|\mbox{\boldmath $A$}(\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i)\|/\|\mbox{\boldmath $b$}\|$",'fontsize', 16,'Interpreter','latex')
colorbar, colormap('cool')

saveas(gcf,'../figures/grid_conds_res.fig')

%% plot
[ress_plt,conds_plt] = meshgrid(ress',conds);
figure
subplot(121)
datamat = forerr;
contourf(log10(ress_plt/eps),log10(conds_plt),log10(datamat), 9), clim([-16, 10])
set(gca,'fontsize', 14)
xlabel('$\|\boldmath $b$ - \mbox{\boldmath $A$}\mbox{\boldmath $x$}\|$', 'Interpreter','latex', 'fontsize', 14)
ylabel('$\kappa$', 'Interpreter','latex', 'fontsize', 14)
title("Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$",'fontsize', 14,'Interpreter','latex')
colorbar, colormap('cool')

subplot(122)
datamat = bwerr;
contourf(log10(ress_plt/eps),log10(conds_plt),log10(datamat), 5),
set(gca,'fontsize', 14)
xlabel('$\|b-Ax\|/u$', 'Interpreter','latex', 'fontsize', 14)
ylabel('$\kappa$', 'Interpreter','latex', 'fontsize', 14)
title("Backward error $\mbox{BE}_\infty(\mbox{\boldmath $\widehat{x}$}_i)$",'fontsize', 14,'Interpreter','latex')
colorbar, colormap('cool')

%% plot
[ress_plt,conds_plt] = meshgrid(ress',conds);
figure('Position', [353   452   774   303])

subplot(121)
datamat = itcountCR;
min_data = 0;
log_cutoff= 40;
max_log_scale =log10(120/40); %number of orders of magnitude spanned by the log scale
data_to_plot = (datamat<log_cutoff).*(datamat-min_data)./log_cutoff + ...
    (datamat>log_cutoff).*(datamat<log_cutoff*10^max_log_scale).*((log10(datamat)-log10(log_cutoff))/max_log_scale*0.2 + 1)+...
    (datamat>log_cutoff*10^max_log_scale)*1.2  + 1;
contourf(log10(ress_plt),log10(conds_plt),data_to_plot,15)
set(gca,'fontsize', 14)
xlabel('$\|\mbox{\boldmath $b$} - \mbox{\boldmath $A$}\mbox{\boldmath $x$}\|$', 'Interpreter','latex', 'fontsize', 14)
ylabel('$\kappa$', 'Interpreter','latex', 'fontsize', 14)
title("Iterations ($m=4000,n=50$)",'fontsize', 14,'Interpreter','latex')

colorbar(gca, 'Xtick', [1/4,2/4,3/4,1,1.2^0.5,1.2] + 1,'Xticklabel',{'10','20','30','40','80','120'}) 
clim([1 2.2])
colormap([cool(50); [ones(9,1),linspace(0.1,1,9)'.^0.5,ones(9,1)]])



subplot(122)
[nvec_plt,mvec_plt] = meshgrid(nvec',mvec);
datamat = itcountMN;
datamat(datamat == 0) = NaN;

min_data = 20;
log_cutoff= 30;
max_overall =80;
tic_distance_linear =10;
max_log_scale =log10(max_overall/log_cutoff); %number of orders of magnitude spanned by the log scale

data_to_plot = (datamat<=log_cutoff).*(datamat-min_data)./log_cutoff + ...
    (datamat>log_cutoff).*(datamat<log_cutoff*10^max_log_scale).*((log10(datamat)-log10(log_cutoff))/max_log_scale*0.2 + ((datamat-min_data)./log_cutoff));

contourf(log10(nvec_plt),log10(mvec_plt),data_to_plot,15)
set(gca,'fontsize', 14)
xlabel('$n$', 'Interpreter','latex', 'fontsize', 14)
ylabel('$m$', 'Interpreter','latex', 'fontsize', 14)
title("Iterations ($\kappa=10^8, \|b-Ax\|/u \approx 10^{12}$)",'fontsize', 14,'Interpreter','latex')
xlim([log10(nvec(1)) log10(nvec(end))])
%colorbar(gca, 'Xtick', [(tic_distance_linear:tic_distance_linear:log_cutoff)/log_cutoff,1.2],...
%    'Xticklabel',[cellstr(num2str((tic_distance_linear:tic_distance_linear:log_cutoff)'))', ...
%    {num2str(max_overall)}])
colorbar(gca, 'Xtick', [0 0.1667 1/3  0.7253 1.1042 1.4747  1.8394 2.2],...
    'Xticklabel',{'20', '25', '30', '40', '50', '60', '70', '80'})
clim([0 2.2])
colormap([cool(50); [ones(9,1),linspace(0.1,1,9)'.^0.5,ones(9,1)]])

