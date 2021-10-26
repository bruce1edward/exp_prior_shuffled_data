rng('default')
addpath('..\..\functions')
k_n = [0.2 0.25 0.3 0.35 0.4 0.45 0.5];nmethod = 6;
iters = 100;iter = 400; % number of iteration for EM alogrithm
mcmc_steps = 4000; burn_steps = 2000; 
error_est = zeros(numel(k_n),iters,nmethod);
r_square = zeros(numel(k_n),iters,nmethod);
dev = @(Y,X,beta) 2*nansum(Y.*log(Y./exp(X*beta)) - (Y - exp(X*beta)));
diag_LL1 = zeros(numel(k_n),iters);
diag_C1 = zeros(numel(k_n),iters);
for i = 1:numel(k_n)
for j = 1:iters
n = 1000; d= 20; r = 50; k = floor(r*k_n(i)); b = 3; 
[X,Y_P,beta,Pi,blockindex,Q] = generate_distribution_block_GLM(n, d, r, k, b);
%Naive 
glm = fitglm(X,Y_P,'linear','Distribution','poisson'); 
beta_naive = glm.Coefficients.Estimate;
%Oracle
glm1 = fitglm(X(Pi,:),Y_P,'linear','Distribution','poisson'); 
beta_oracle = glm1.Coefficients.Estimate;
%EM
beta_EM = EM_mal_Block_GLM2(Y_P, X, iter, mcmc_steps, burn_steps, beta_naive, blockindex, 0);
%EMM
theta = Choose_theta(r, k);
beta_EMM = EM_mal_Block_GLM2(Y_P, X, iter, mcmc_steps, burn_steps, beta_naive, blockindex, theta);
%LL
[beta_LL, norm_score_LL, converged_LL] = LL_Armijo(n,d,X,Q,Y_P,beta_naive);
%C
[beta_C, norm_score_C, converged_C] = Chamber_Armijo(n,d,X,Q,Y_P,beta_naive);
diag_LL(i,j) = converged_LL;
diag_C(i,j) = converged_C;

error_est(i,j,1) = norm(beta_naive - beta)/b;
error_est(i,j,2) = norm(beta_oracle  - beta)/b;
error_est(i,j,3) = norm(beta_EM  - beta)/b;
error_est(i,j,4) = norm(beta_EMM  - beta)/b;
error_est(i,j,5) = norm(beta_LL  - beta)/b;
error_est(i,j,6) = norm(beta_C  - beta)/b;
end
end

save('Block_GLM.mat');


