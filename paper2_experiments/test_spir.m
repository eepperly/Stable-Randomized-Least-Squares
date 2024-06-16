addpath('../utils')
addpath('../code')

rng(3423390)
m = 4000;
n = 50;
d = 12*n;
[A,b,x,r] = random_ls_problem(4000,50,1e12,1e-6);
colors
trials = 25;
lastiter = 25;

real_run = true;
if real_run
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);backward_error_ls(A,b,y)/norm(A,'fro')];
else
    [~,S,V] = svd(A,"econ"); %#ok<*UNRCH>
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);kw_estimate(A,b,y,Inf,S,V)/norm(A,'fro')];
end

lsqr_runs = nan * zeros(25,3,trials);
cg_runs  = nan * zeros(25,3,trials);
foss_runs  = nan * zeros(25,3,trials);
for trial = 1:trials
    [~,spir_lsqr] = spir(A,b,d,'halfadaptive',summary,'lsqr',true);
    [~,spir_cg] = spir(A,b,d,'halfadaptive',summary,'cg',true);
    [~,foss] = fossils(A,b,d,'halfadaptive',summary,true);
    lsqr_runs(1:size(spir_lsqr,1),1:3,trial) = spir_lsqr;
    cg_runs(1:size(spir_cg,1),1:3,trial) = spir_cg;
    foss_runs(1:size(foss,1),1:3,trial) = foss;
end

spir_lsqr = [];
spir_cg_upper = [];
spir_lsqr_upper = [];
spir_cg = [];
spir_lsqr_lower = [];
spir_cg_lower = [];
foss = [];
foss_upper = [];
foss_lower = [];
for i = 1:(lastiter+1)
    lsqrstuff = lsqr_runs(i,:,:);
    badinds = [];
    for trial = 1:trials
        if isnan(lsqrstuff(1,1,trial))
            badinds = [badinds trial]; %#ok<*AGROW>
        end
    end
    lsqrstuff(:,:,badinds) = [];
    if ~isempty(lsqrstuff)
        spir_lsqr(i,:) = median(lsqrstuff,3); %#ok<*SAGROW>
        spir_lsqr_upper(i,:) = quantile(lsqrstuff,0.8,3);
        spir_lsqr_lower(i,:) = quantile(lsqrstuff,0.2,3);
    end

    cgstuff = cg_runs(i,:,:);
    badinds = [];
    for trial = 1:trials
        if isnan(cgstuff(1,1,trial))
            badinds = [badinds trial]; 
        end
    end
    cgstuff(:,:,badinds) = [];
    if ~isempty(cgstuff)
        spir_cg(i,:) = median(cgstuff,3);
        spir_cg_upper(i,:) = quantile(cgstuff,0.8,3);
        spir_cg_lower(i,:) = quantile(cgstuff,0.2,3);
    end

    fossstuff = foss_runs(i,:,:);
    badinds = [];
    for trial = 1:trials
        if isnan(fossstuff(1,1,trial))
            badinds = [badinds trial];
        end
    end
    fossstuff(:,:,badinds) = [];
    if ~isempty(fossstuff)
        foss(i,:) = median(fossstuff,3); 
        foss_upper(i,:) = quantile(fossstuff,0.8,3);
        foss_lower(i,:) = quantile(fossstuff,0.2,3);
    end
end

[Q,R] = qr(A,'econ');
y = R\(Q'*b);
qr_vals = summary(y);

%% Plot

close all
for j = 1:3
    figure(j)
    plot_shaded(0:(size(spir_lsqr,1)-1), spir_lsqr(:,j), ...
        spir_lsqr_lower(:,j), spir_lsqr_upper(:,j), purple,...
        'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', purple, ...
        'Marker', 's');
    hold on
    plot_shaded(0:(size(spir_cg,1)-1), spir_cg(:,j), ...
        spir_cg_lower(:,j), spir_cg_upper(:,j), pink,...
        'LineWidth', 1, 'MarkerSize', 10,...
        'MarkerFaceColor', pink, 'Marker', 'v', 'LineStyle', ':');
    plot_shaded(0:(size(foss,1)-1), foss(:,j), ...
        foss_lower(:,j), foss_upper(:,j), orange,...
        'LineWidth', 1, 'MarkerSize', 10,...
        'MarkerFaceColor', orange, 'Marker', '^', 'LineStyle', '--'); 
    set(gca, 'YScale', 'log')
    yline(qr_vals(j),'k:', 'LineWidth', 3)
    xlabel('Iteration $i$')
    if j == 1
        ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
        legend('SPIR (LSQR)','','SPIR (CG)','','FOSSILS','','QR')
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