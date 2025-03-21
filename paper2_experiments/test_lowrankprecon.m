%% Set up problem

addpath('../code')
addpath('../utils')
colors

rng(32439)

m = 4000;
n = 50;

trials = 25; 

cond_A = 10^20;
res_size = 10^-12;

[A,b,x,r,S,V] = random_ls_problem(m,n,cond_A,res_size);

real_run = false;
if real_run
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);backward_error_ls(A,b,y)/norm(A,'fro')];
else
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);kw_estimate(A,b,y,Inf,S,V)/norm(A,'fro')];
end

%% FOSSILS
'Fossils 1'
[~,fossil]=fossils(A,b,12*n,[],summary,true);

%% SPIR
'SPIR 1'
[~,spirs]=spir(A,b,12*n,[],summary,'lsqr',true);

%% FOSSILS
'FOSSILS 1'
[~,fossil_lowrank]=fossils(A,b,12*n,[],summary,true, false, true);

%% SPIR
'SPIR 1'
[~,spirs_lowrank]=spir(A,b,12*n,[],summary,'lsqr',true, false, true);

%% QR

[Q,R] = qr(A,'econ');
y = R\(Q'*b);
qr_vals = summary(y);

%% Plot

figure
for j = 1:3
    subplot(1,3,j)
    semilogy(0:(length(fossil(:,j))-1), fossil(:,j), '^-', 'LineWidth', 1, 'Color',...
        orange, 'MarkerSize', 10, 'MarkerFaceColor', orange); 
    hold on
    semilogy(0:(length(spirs(:,j))-1), spirs(:,j), 'v-', 'LineWidth', 1, 'Color',...
        pink, 'MarkerSize', 10, 'MarkerFaceColor', pink); 
    semilogy(0:(length(fossil_lowrank(:,j))-1), fossil_lowrank(:,j), 's:', 'LineWidth', 1, 'Color',...
        blue, 'MarkerSize', 10, 'MarkerFaceColor', blue); 
    semilogy(0:(length(spirs_lowrank(:,j))-1), spirs_lowrank(:,j), '*:', 'LineWidth', 1, 'Color',...
        purple, 'MarkerSize', 10, 'MarkerFaceColor', purple); 
    yline(qr_vals(j),'k:', 'LineWidth', 3) 
    hold off
    xlabel('Iteration $i$')
    if j == 1
        ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
        
	if real_run
            saveas(gcf,'../figs/compare_forward.fig')
            saveas(gcf,'../figs/compare_forward.png')
	end
    elseif j == 2
        ylabel('Residual error $\|\mbox{\boldmath $A$}(\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i)\|/\|\mbox{\boldmath $b$}\|$')
	if real_run
            saveas(gcf,'../figs/compare_residual.fig')
            saveas(gcf,'../figs/compare_residual.png')
	end
    elseif j == 3
        ylabel('Backward error $\mbox{BE}_\infty(\mbox{\boldmath $\widehat{x}$}_i)$')
	if real_run
            saveas(gcf,'../figs/compare_backward.fig')
            saveas(gcf,'../figs/compare_backward.png')
	end
    end
    legend('FOSSILS (explicit reg.)',...
            'SPIR (explicit reg.)',...
            'FOSSILS (low-rank precon.)',...
            'SPIR (low-rank precon.)',...
            'QR')
    set(gca,'FontSize',14)
end
sgtitle(['$A$ is ', num2str(m), '$\times$', num2str(n), ' with $\kappa_2(A) = $', num2str(cond(A), '%.1e'), ' and $\|r^*\| = $', num2str(res_size)])

%% Save

if real_run
    save('../data/results_sketch_poly.mat', 'lsqrwarm', 'lsqrcold', 'itsk', 'fossil', 'mom', 'qr_vals', 'trials')
end
