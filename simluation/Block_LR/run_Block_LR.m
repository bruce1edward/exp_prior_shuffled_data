rng('default')
addpath('..\..\functions')
k_n = [0.2 0.25 0.3 0.35 0.4 0.45 0.5];nmethod = 6;
iters = 100;iter = 400; % number of iteration for EM alogrithm
mcmc_steps = 4000; burn_steps = 2000; 
error_est = zeros(numel(k_n),iters,nmethod);
for i = 1:numel(k_n)
for j = 1:iters
n = 1000; d= 20; r = 50; K = floor(r*k_n(i)); sigma = 1; b = 3; 
[X,Y,Y_P,beta,Pi,blockindex,Q] = generate_distribution_block(n, d, r, K, sigma, b);
%Naive 
beta_naive = X\Y_P;
%Oracle
beta_oracle = X\Y;
%LL
QX = Q*X;
beta_LL = QX\Y_P;
%Chamber
beta_C = (X'*Q*X)\(X'*Y_P);
%EM
[beta_EM, sigma_EM] = EM_mal_Block(Y_P, X, iter, mcmc_steps, burn_steps, 0, beta_naive, blockindex);
%EMM
theta = Choose_theta(r, K);
[beta_EMM, sigma_EMM] =  EM_mal_Block(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, blockindex);

error_est(i,j,1) = norm(beta_naive - beta)/b;
error_est(i,j,2) = norm(beta_oracle  - beta)/b;
error_est(i,j,3) = norm(beta_EM  - beta)/b;
error_est(i,j,4) = norm(beta_EMM  - beta)/b;
error_est(i,j,5) = norm(beta_LL  - beta)/b;
error_est(i,j,6) = norm(beta_C  - beta)/b;
end
end
save('Block_LR.mat');

