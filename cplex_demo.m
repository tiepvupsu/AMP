function sparse_recovery_demo_pos1(Atoms, len, TrueSparsityLevel, pos)
    % reweighted_l1, Elastic net, sparsa, icr, amp 
clc
% clear all
close all
myinit();
addpath('utils');
addpath('ICR');
addpath('Enet');
addpath('AMP');
addpath('NNOMP');
ind=1;
if nargin == 0 
%     Atoms = 64; len = 32; TrueSparsityLevel = 10;
    Atoms = 512; len = 256; TrueSparsityLevel = 100;
end 
if nargin < 4
    pos = 1;
end 
ntest = 10;
fn = fullfile('results', strcat('SR_d', num2str(len), '_k', num2str(Atoms), ...
        '_L', num2str(TrueSparsityLevel), '_pos_', num2str(pos), '.mat'));
disp(fn);

opts.pos = pos; % positive constraint in the solution 

ICR    = zeros(ntest, 5);
ENET   = zeros(ntest, 5);
RWL1   = zeros(ntest, 5);
SPARSA = zeros(ntest, 5);
AMP    = zeros(ntest, 5);
OMP    = zeros(ntest, 5);

%% evaluate the result: 
%           MSE, cost_fn (spike and slab), sparsity level, support match 
function res = evaluate_res(x1, t1)
    x0 = xg;
    MSE = mse(x1, x0);
    cost_fn = cost_spike_slab(y, A, x1, lambda, Rho);
    SL = sum(x1 ~= 0);
    sm = SM(x1, x0);
    res = [t1, MSE, cost_fn, SL, sm];
end 

opts.display = 0;
for i = 1: ntest
    i
    fn = strcat('cplex_', num2str(i), '.mat');
    load(fn);
    
    if i == 10
    A = A_cell{i};
    Rho = Rho_cell{i};
    y = y_cell{i};
    xg = xg_cell{i};
    lambda = lambda_cell{i};
    else 
        A = A;
        Rho = Rho;
        y = y;
        xg = xg;
        lambda = lambda;
    end 
    Lambda =0.0002; 
    Sigma = 0.018;
    Kappa = 0.47 * ones(Atoms,1);
%     lambda = Lambda;
    
    x0 = xg;
    
    %% ========================== reweighted l1 =============================
    tic
    x1 = reweighted_l1(y,A, Sigma, opts);
    t1 = toc;
    
    %% ========================= Elastic net=========================    
    tic;
    x2 = Enet(y, A, lambda, 0.1, opts);
    t2 = toc;

    %% ===================== Sparsa =============================    
    tic;
    x3 = my_sparsa_sc(y, A, lambda, Rho, opts);
    t3 = toc;
    %% ========================== ICR =============================
    AtA = A'*A;
    y1 =                                                                                                                                                                                    (y);
    tic;    
%         [x4,~,~,~] = ICR_Func(y1,A,AtA,'lambda',Lambda,...
%             'rho',Rho,'GroundTruth',x0, 'ALGORITHM', 1-pos, 'VERBOSE', 0 );
    [x4,~,~,~] = ICR_Func_modified(y,A,AtA,'lambda',Lambda,...
        'rho',Rho,'GroundTruth',x0, 'ALGORITHM', 1-pos, 'VERBOSE', 0 );
%     x4 = rand(size(x0));
    t4 = toc;

    %% ===================== AMP ==============================
    tic 
    x5 = AMP_chol(y, A, lambda, Rho, opts);
    t5 = toc;
    % ---------------------- end of AMP --------------------------
    % 
    %% ===================== NNOMP ==============================
    param.L=10; % not more than 10 non-zeros coefficients
    param.eps=0.01; % squared norm of the residual should be less than 0.1
    param.numThreads=-1; % number of processors/cores to use; the default choice is -1
                        % and uses all the cores of the machine
    param.pos = pos;
    tic
%     x6=full(mexOMP(y,A,param));
%     x6 = NNOMP(y, A, 10);
    opts.printEvery = 100;
    x6 = CoSaMP( A, y, 10, [], opts )';
    % x1 = reweighted_l1(y,A, Sigma, opts);
    % t=toc
    t6 = toc;
    x6 = full(x6');
    %------------------    
    res1 = evaluate_res(x1, t1);
    res2 = evaluate_res(x2, t2);
    res3 = evaluate_res(x3, t3);
    res4 = evaluate_res(x4, t4);
    res5 = evaluate_res(x5, t5);
    res6 = evaluate_res(x6, t6);

    fprintf('                | Time    |   MSE   |  cost  |   SL   |   SM  |\n');
    % fprintf('')
    fprintf('Reweighted l1:')
    disp(res1);
    fprintf('Elastic Net  :');
    disp(res2);
    fprintf('SPARSA       :');
    disp(res3);
    fprintf('ICR          :');
    disp(res4);
    fprintf('CoSaMP        :');
    disp(res6);
    fprintf('AMP          :');
    disp(res5);
    
    fprintf('cost0 = %f | SL = %d\n', cost_spike_slab(y, A, x0, lambda, Rho), ...
                                        sum(x0 ~= 0));
%     imagesc([x0, x1, x3]);
%     pause
    RWL1(i,:) = res1;
    ENET(i,:) = res2;
    SPARSA(i,:) = res3;
    ICR(i,:) = res4;
    AMP(i,:) = res5;
    OMP(i,:) = res6;
end

fprintf('\n\nOverall:\n')
fprintf('                | Time    |   MSE   |  cost  |   SL   |   SM  |\n');
fprintf('Reweighted l1:');
disp(mean(RWL1));
fprintf('Elastic Net  :');
disp(mean(ENET));
fprintf('SPARSA       :');
disp(mean(SPARSA));
fprintf('ICR          :');
disp(mean(ICR));
fprintf('CoSaMP        :');
disp(mean(OMP));
fprintf('AMP          :');
disp(mean(AMP));


rwl1 = mean(RWL1);
enet = mean(ENET);
sparsa = mean(SPARSA);
icr = mean(ICR);
amp = mean(AMP);
omp = mean(OMP);
% cosamp = rwl1;
fn = 'cplex_res.mat';
save(fn, 'rwl1', 'enet', 'sparsa', 'icr', 'amp', 'omp');

end 


function res = sm2(x1, x2)
    s1 = find(x1);
    s2 = find(x2);
    res = numel(intersect(s1, s2)) /numel(union(s1, s2));
end 
