%% Set up problem

addpath('../code')
addpath('../utils')
colors

rng(235029)

m = 4000;
n = 50;

difficulties = logspace(0,16,21);
real_run = true;

if real_run
    be = @(y,A,b,S,V) backward_error_ls(A,b,y)/norm(A,'fro');
else
    be = @(y,A,b,S,V) kw_estimate(A,b,y,Inf,S,V)/norm(A,'fro');
end

spre = zeros(length(difficulties),3);
imom = zeros(length(difficulties),3);
foss = zeros(length(difficulties),3);
foss1 = zeros(length(difficulties),3);
spirs = zeros(length(difficulties),3);
qrs  = zeros(length(difficulties),3);

for d_idx = 1:length(difficulties)
    difficulty = difficulties(d_idx)
    [A,b,x,r,S,V] = random_ls_problem(m,n,difficulty,eps*difficulty);

    % spre 
    xhat = sketch_and_precondition(A,b,12*n);
    spre(d_idx,1) = norm(x-xhat);
    spre(d_idx,2) = norm(A*(x-xhat));
    spre(d_idx,3) = be(xhat,A,b,S,V);

    % spir
    xhat = spir(A,b,12*n);
    spirs(d_idx,1) = norm(x-xhat);
    spirs(d_idx,2) = norm(A*(x-xhat));
    spirs(d_idx,3) = be(xhat,A,b,S,V);

    % imom
    xhat = iterative_sketching(A,b,[],[],[],[],'optimal','optimal');
    imom(d_idx,1) = norm(x-xhat);
    imom(d_idx,2) = norm(A*(x-xhat));
    imom(d_idx,3) = be(xhat,A,b,S,V);

    % foss
    xhat = fossils(A,b);
    foss(d_idx,1) = norm(x-xhat);
    foss(d_idx,2) = norm(A*(x-xhat));
    foss(d_idx,3) = be(xhat,A,b,S,V);

    % foss1
    xhat = fossils(A,b,[],100);
    foss1(d_idx,1) = norm(x-xhat);
    foss1(d_idx,2) = norm(A*(x-xhat));
    foss1(d_idx,3) = be(xhat,A,b,S,V);

    % qrs
    [Q,R] = qr(A,0); xhat = R\(Q'*b);
    qrs(d_idx,1) = norm(x-xhat);
    qrs(d_idx,2) = norm(A*(x-xhat));
    qrs(d_idx,3) = be(xhat,A,b,S,V);
end

%% Plot

close all
for j = 1:3
    figure(j)
    loglog(difficulties, spre(:,j), 's-', 'LineWidth', 1, 'Color',...
        purple, 'MarkerSize', 10, 'MarkerFaceColor', purple); hold on
    loglog(difficulties, imom(:,j), 'o-', 'LineWidth', 1, 'Color',...
        blue,'MarkerSize', 10, 'MarkerFaceColor', blue); 
    loglog(difficulties, foss1(:,j), '-.', 'LineWidth', 4, 'Color', yellow); 
    loglog(difficulties, foss(:,j), '^-', 'LineWidth', 1, 'Color',...
        orange, 'MarkerSize', 10, 'MarkerFaceColor', orange); 
    loglog(difficulties, spirs(:,j), 'v:', 'LineWidth', 1, 'Color',...
        pink, 'MarkerSize', 10, 'MarkerFaceColor', pink); 
    loglog(difficulties, qrs(:,j), '--', 'LineWidth', 2, 'Color', black); 
    axis([1e0 1e16 -Inf Inf])
    xlabel('Difficulty = $\kappa = \|\mbox{\boldmath $b$}-\mbox{\boldmath $Ax$}\|/u$')
    if j == 1
        ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}\|/\|\mbox{\boldmath $x$}\|$')
        legend({'Sketch\&Pre','Iter Sketch + Mom','FOSSILS (1 step)','FOSSILS','SPIR','QR'},'Location','northwest')
	if real_run
            saveas(gcf,'../figs/compare_2_forward.fig')
            saveas(gcf,'../figs/compare_2_forward.png')
	end
    elseif j == 2
        ylabel('Residual error $\|\mbox{\boldmath $A$}(\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$})\|/\|\mbox{\boldmath $b$}\|$')
	if real_run
            saveas(gcf,'../figs/compare_2_residual.fig')
            saveas(gcf,'../figs/compare_2_residual.png')
	end
    elseif j == 3
        ylabel('Backward error $\mbox{BE}_\infty(\mbox{\boldmath $\widehat{x}$})$')
	if real_run
            saveas(gcf,'../figs/compare_2_backward.fig')
            saveas(gcf,'../figs/compare_2_backward.png')
	end
    end
end

%% Save

if real_run
    save('../data/results_compare_2.mat', 'spre', 'spirs', 'itsk', 'foss', 'foss1', 'imom', 'qrs')
end
