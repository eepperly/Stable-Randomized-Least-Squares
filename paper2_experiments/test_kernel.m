%% Load SUSY data

addpath('../code')
addpath('../utils')

load('../data/susy.mat')
X = X(1:1e6,:); b = b(1:1e6);
bandwidth = 4;
S = randsample(1e6,1e3);
A = exp(-pdist2(X,X(S,:),"euclidean").^2 / (2*bandwidth^2));

sizes = round(logspace(1,3,11));

real_run = true;

%% Run trials

dir = zeros(length(sizes), 1);
moms = zeros(length(sizes), 1);
foss = zeros(length(sizes), 1);

for i = 1:11
    sizes(i)
    AA = A(:,1:sizes(i));

    % QR
    tic; x_qr = AA\b; dir(i) = toc;
    dir(i)

    % Iterative sketching with momentum
    tic; x = iterative_sketching(AA,b,12*sizes(i),[],[],[],'optimal','optimal');
    moms(i) = toc;
    moms(i)

    % FOSSILS
    tic; x = fossils(AA,b);
    foss(i) = toc;
    foss(i)
end

if real_run
    save('../data/results_timing.mat', 'dir', 'moms', 'foss')
end

%% Plot
addpath('../utils/')
colors

close all
figure(1)
loglog(sizes, dir(:,1), '--', 'LineWidth', 2, 'Color', black); hold on
loglog(sizes, moms(:,1), 'o-', 'LineWidth', 1, 'Color',...
        blue,'MarkerSize', 10, 'MarkerFaceColor', blue); 
loglog(sizes, foss(:,1), '^-', 'LineWidth', 1, 'Color',...
        orange, 'MarkerSize', 10, 'MarkerFaceColor', orange); 

xlabel('Number of columns $n$')
ylabel('Time (sec)')
% axis([-Inf Inf 6e-2 2e2])
legend({'\texttt{mldivide}','Iter Sketch + Mom','FOSSILS',},'Location','northwest')
saveas(gcf,'../figs/kernel.fig')
saveas(gcf,'../figs/kernel.png')