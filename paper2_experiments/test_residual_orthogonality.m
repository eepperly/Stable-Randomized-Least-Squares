m = 4000;
n = 50;
[A,b] = random_ls_problem(m,n,1e12,1e-3);

trials = 100;
qrs = zeros(trials,1);
imom = zeros(trials,1);
spre = zeros(trials,1);
foss = zeros(trials,1);

for trial = 1:trials
    % QR
    [Q,R] = qr(A,"econ");
    x = R\(Q'*b);
    qrs(trial) = norm(A'*(b - A*x));
    
    % Iterative sketching with momentum
    x = iterative_sketching(A,b,[],[],[],[],"optimal","optimal");
    imom(trial) = norm(A'*(b - A*x));
    
    % Sketch-and-precondition
    x = sketch_and_precondition(A,b);
    spre(trial) = norm(A'*(b - A*x));
    
    % FOSSILS
    x = fossils(A,b);
    foss(trial) = norm(A'*(b - A*x));
end

disp('QR'); median(qrs)
disp('imom'); median(imom)
disp('spre'); median(spre)
disp('FOSSILS'); median(foss)