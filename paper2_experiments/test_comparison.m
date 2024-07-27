%% Set up problem

addpath('../code')
addpath('../utils')
colors

rng(32439)

m = 4000;
n = 50;

trials = 25; 

cond_A = 10^12;
res_size = 10^-6;

[A,b,x,r,S,V] = random_ls_problem(m,n,cond_A,res_size);

real_run = true;
if real_run
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);backward_error_ls(A,b,y)/norm(A,'fro')];
else
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);kw_estimate(A,b,y,Inf,S,V)/norm(A,'fro')];
end

%% LSQR Warm

[~,lsqrwarm]=sketch_and_precondition(A,b,12*n,trials,summary,true,'lsqrwarm',true);

%% LSQR Cold

[~,lsqrcold]=sketch_and_precondition(A,b,12*n,trials,summary,true,'lsqrcold',true);

%% Iterative Sketching with Momentum

[~,mom]=iterative_sketching(A,b,12*n,trials,summary,true,'optimal','optimal',true);

%% FOSSILS

[~,fossil]=fossils(A,b,12*n,[],summary,true);

%% SPIR

[~,spirs]=spir(A,b,12*n,[],summary,[],true);

%% QR

[Q,R] = qr(A,'econ');
y = R\(Q'*b);
qr_vals = summary(y);

%% Plot

close all
for j = 1:3
    figure(j)
    semilogy(0:trials, lsqrcold(:,j), '-.', 'LineWidth', 4, 'Color', yellow); hold on
    semilogy(0:trials, mom(:,j), 'o-', 'LineWidth', 1, 'Color',...
        blue,'MarkerSize', 10, 'MarkerFaceColor', blue); 
    semilogy(0:trials, lsqrwarm(:,j), 's-', 'LineWidth', 1, 'Color',...
        purple, 'MarkerSize', 10, 'MarkerFaceColor', purple);
    semilogy(0:(length(fossil(:,j))-1), fossil(:,j), '^-', 'LineWidth', 1, 'Color',...
        orange, 'MarkerSize', 10, 'MarkerFaceColor', orange); 
    semilogy(0:(length(spirs(:,j))-1), spirs(:,j), 'v:', 'LineWidth', 1, 'Color',...
        pink, 'MarkerSize', 10, 'MarkerFaceColor', pink); 
    yline(qr_vals(j),'k:', 'LineWidth', 3)
    xlabel('Iteration $i$')
    if j == 1
        ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
        legend('Sketch\&Pre ($\mbox{\boldmath $x$}_0=\mbox{\boldmath $0$}$)',...
            'Iter Sketch + Mom',...
            'Sketch\&Pre',...
            'FOSSILS',...
            'SPIR',...
            'QR')
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
end

%% Save

if real_run
    save('../data/results_sketch_poly.mat', 'lsqrwarm', 'lsqrcold', 'itsk', 'fossil', 'mom', 'qr_vals', 'trials')
end
