rng('default')
addpath('..\..\functions')
k_n = [0.2 0.25 0.3 0.35 0.4 0.45 0.5];nmethod = 5;
iters = 100;iter = 400; % number of iteration for EM alogrithm
mcmc_steps = 8000; burn_steps = 4000; 
error_est = zeros(numel(k_n),iters,nmethod);
r_square = zeros(numel(k_n),iters,nmethod);
dev = @(Y,X,beta) 2*nansum(Y.*log(Y./exp(X*beta)) - (Y - exp(X*beta))); 

for i = 1:numel(k_n)
for j = 1:iters
n = 1000; d= 20; K = floor(n*k_n(i)); b = 1;
[X,Y_P,beta,Pi,inv_Pi] = generate_distribution_sparse_GLM(n, d, K, b);
%Naive 
glm = fitglm(X,Y_P,'linear','Distribution','poisson'); 
beta_naive = glm.Coefficients.Estimate;
%Oracle
glm1 = fitglm(X(Pi,:),Y_P,'linear','Distribution','poisson'); 
beta_oracle = glm1.Coefficients.Estimate;
%EM
order = 1:n;
beta_EM = EM_mal_GLM(n, Y_P, X, iter, mcmc_steps, burn_steps, 0, beta_naive, order);
%EMM
order = 1:n;
theta = Choose_theta(n,K);
beta_EMM = EM_mal_GLM(n, Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);
%robust
lambda = 0.5*sqrt(mean(Y_P))*sqrt(log(n+d)/n);
[beta_robust, xi_robust] = AEA_Poi(X,Y_P,lambda,beta_naive);

error_est(i,j,1) = norm(beta_naive - beta)/b;
error_est(i,j,2) = norm(beta_oracle  - beta)/b;
error_est(i,j,3) = norm(beta_EM  - beta)/b;
error_est(i,j,4) = norm(beta_EMM  - beta)/b;
error_est(i,j,5) = norm(beta_robust  - beta)/b;
end
end

save('Hamming_GLM.mat');

