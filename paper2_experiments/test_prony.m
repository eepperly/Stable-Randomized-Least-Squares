addpath('../code/')
H = tfim(16,1);
dt = 0.01;
psi = zeros(size(H,1),1); psi(1) = 1; last_psi = psi;
batch_size = 100;
num_batches = 7.6e3;

real_run = true;
est = @(A,b,x) kw_estimate(A,b,x,Inf,"sketched") / norm(A,'fro');

try load('../data/prony.mat'); catch; end
if ~exist('inprods','var') || length(inprods) < batch_size*num_batches
    fprintf('Generating data...\n')
    inprods = zeros(batch_size*num_batches,1);
    hadtest = @(x) 2*(rand(length(x),1) < (x+1)/2)-1;
    for batch = 1:num_batches
        disp(batch)
        evolves = expmv(1i*H,last_psi,dt*(1:batch_size));
        inprods(((batch-1)*batch_size+1):(batch*batch_size)) = psi'*evolves;
        last_psi = evolves(:,end);
        if mod(batch, 10) == 0
            save('../data/prony.mat', 'inprods')
        end
    end
end

measurements = inprods + 1e-6/sqrt(2) * (randn(size(inprods))...
    + 1i*randn(size(inprods)));

m = 7.5e5;
ns = round(logspace(log10(50),3,11));
foss = zeros(length(ns),1);
spirs = zeros(length(ns),1);
itsk = zeros(length(ns),1);
dir = zeros(length(ns),1);

fprintf('n\tfoss\t\tspir\t\tism\t\tdir\n')
for n_idx = 1:length(ns)
    n = ns(n_idx);
    A = toeplitz(measurements(n:(n+m-1)), measurements(n:-1:1).');
    b = measurements((n+1):(n+m));

    foss(n_idx) = timeit(@() fossils(A,b));
    spirs(n_idx) = timeit(@() spir(A,b));
    itsk(n_idx) = timeit(@() iterative_sketching(A,b,12*n,[],[],[],...
        'optimal','optimal'));
    dir(n_idx) = timeit(@() A\b);

    fprintf('%d\t%e\t%e\t%e\t%e\n', n, foss(n_idx), spirs(n_idx),...
        itsk(n_idx), dir(n_idx))
    if real_run
        save('../data/prony.mat','foss','itsk','dir','ns','m')
    end
end

%% Plot
addpath('../utils/')
colors

close all
figure(1)
loglog(ns, dir(:,1), '--', 'LineWidth', 2, 'Color', black); hold on
loglog(ns, itsk(:,1), 'o-', 'LineWidth', 1, 'Color',...
        blue,'MarkerSize', 10, 'MarkerFaceColor', blue); 
loglog(ns, foss(:,1), '^-', 'LineWidth', 1, 'Color',...
        orange, 'MarkerSize', 10, 'MarkerFaceColor', orange); 
loglog(ns, spirs(:,1), 'v:', 'LineWidth', 1, 'Color',...
        pink, 'MarkerSize', 10, 'MarkerFaceColor', pink); 

xlabel('Number of columns $n$')
ylabel('Time (sec)')
axis([-Inf Inf 4e-2 2e2])
legend({'\texttt{mldivide}','Iter Sketch + Mom','FOSSILS','SPIR'},'Location','northwest')
if real_run
    saveas(gcf,'../figs/prony.fig')
    saveas(gcf,'../figs/prony.png')
end