% All rectangular matrices in suitesparse collection with at least 1e5 rows
% and at most 2e4 columns that are numerically full-rank
addpath('../code/')
fnames = {"bibd_20_10","bibd_22_8","EternityII_A","EternityII_E",...
    "EternityII_Etilde","lp_osa_30","lp_osa_60","nw14","rail2586",...
    "rail4284","spal_004","stat96v1",...
    "mesh_deform","ch7-9-b3","ch8-8-b3","shar_te2-b2"};
idx = [8,1,2,9,10,6,12,3,13,5,11,7,4,16,14,15];
fnames = fnames(idx); % sorted by number of cols

foss = zeros(length(fnames),2);
itsk = zeros(length(fnames),2);
dir  = zeros(length(fnames),2);
ids  = zeros(length(fnames),1);
rows = zeros(length(fnames),1);
cols = zeros(length(fnames),1);
nnzs = zeros(length(fnames),1);
names = cell(length(fnames),1);

for idx = 1:length(fnames)
    try
        load(sprintf("../data/%s.mat",fnames{idx}));
    catch 
        warning("../data/%s.mat not found. Please download this file from the Suitesparse matrix collection. Ignoring and moving on",fnames{idx})
        continue
    end
    ids(idx) = Problem.id;
    A = Problem.A;
    if size(A,1) <= size(A,2)
        A = A';
    end
    size(A,2)
    b = randn(size(A,1),1);

    rows(idx) = size(A,1);
    cols(idx) = size(A,2);
    nnzs(idx) = nnz(A);
    names{idx} = Problem.name;

    tic; x = fossils(A,b,3*size(A,2)); foss(idx,1) = toc;
    foss(idx,2) = kw_estimate(A,b,x,Inf,"sketched")

    tic; x = iterative_sketching(A,b,3*size(A,2),[],[],[],"optimal",...
        "optimal"); itsk(idx,1) = toc;
    itsk(idx,2) = kw_estimate(A,b,x,Inf,"sketched")

    tic; x = A\b; dir(idx,1) = toc;
    dir(idx,2) = kw_estimate(A,b,x,Inf,"sketched")

    fprintf('%d\t%d\t%e\t%e\t%e\n', ids(idx), size(A,1), foss(idx,1),...
        itsk(idx,1), dir(idx,1))
    save('../sparse.mat', 'foss', 'itsk', 'dir')
end