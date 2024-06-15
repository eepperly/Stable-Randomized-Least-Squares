addpath('../utils')
addpath('../code')

m = 4000;
n = 50;
d = 12*n;
[A,b,x,r] = random_ls_problem(4000,50,1e12,1e-6);
trials = 100;
colors

real_run = true;
if real_run
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);backward_error_ls(A,b,y)/norm(A,'fro')];
else
    [~,S,V] = svd(A,"econ");
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);kw_estimate(A,b,y,Inf,S,V)/norm(A,'fro')];
end

[~,spir_lsqr] = spir(A,b,d,[],summary,'lsqr',true);
[~,spir_cg] = spir(A,b,d,[],summary,'cg',true);

[Q,R] = qr(A,'econ');
y = R\(Q'*b);
qr_vals = summary(y);

%% Plot

close all
for j = 1:3
    figure(j)
    semilogy(0:(size(spir_lsqr,1)-1), spir_lsqr(:,j), 's-', 'LineWidth', 1,...
        'Color',purple, 'MarkerSize', 10, 'MarkerFaceColor', purple);
    hold on
    semilogy(0:(size(spir_cg,1)-1), spir_cg(:,j), 'v--', 'LineWidth', 1,...
        'Color', orange, 'MarkerSize', 10, 'MarkerFaceColor', orange); 
    yline(qr_vals(j),'k:', 'LineWidth', 3)
    xlabel('Iteration $i$')
    if j == 1
        ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
        legend('SPIR (LSQR)',...
            'SPIR (CG)',...
            'QR')
	if real_run
            saveas(gcf,'../figs/spir_forward.fig')
            saveas(gcf,'../figs/spir_forward.png')
	end
    elseif j == 2
        ylabel('Residual error $\|\mbox{\boldmath $A$}(\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i)\|/\|\mbox{\boldmath $b$}\|$')
	if real_run
            saveas(gcf,'../figs/spir_residual.fig')
            saveas(gcf,'../figs/spir_residual.png')
	end
    elseif j == 3
        ylabel('Backward error $\mbox{BE}_\infty(\mbox{\boldmath $\widehat{x}$}_i)$')
	if real_run
        saveas(gcf,'../figs/spir_backward.fig')
        saveas(gcf,'../figs/spir_backward.png')
	end
    end
end