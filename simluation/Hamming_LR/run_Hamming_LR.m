rng('default')
addpath('..\..\functions')
k_n = [0.2 0.25 0.3 0.35 0.4 0.45 0.5];nmethod = 6;
iters = 100;iter = 400; % number of iteration for EM alogrithm
mcmc_steps = 8000; burn_steps = 4000; 
error_est = zeros(numel(k_n),iters,nmethod);
r_square = zeros(numel(k_n),iters,nmethod);

for i = 1:numel(k_n)
for j = 1:iters
n = 1000; d= 20; K = floor(n*k_n(i)); sigma = 1; b = 3;
[X,Y,Y_P,beta,Pi,inv_Pi] = generate_distribution_sparse(n, d, K, sigma, b);
X_wo = X(:,2:end);X_centered = X_wo - repmat(mean(X_wo,1), [n 1]);
%Naive 
beta_naive = X\Y_P;
%Oracle
beta_oracle = X\Y;
%Robust
beta_robust = robustfit(X_wo,Y_P,'huber');
%Mixture
Y_P1 = Y_P - mean(Y_P);control = "robust";
[beta_mixture, alpha_mixture, sigma_mixture] = fit_mixture(X_centered, Y_P1, control);
beta_mixture = [mean(Y_P - X_wo*beta_mixture) ; beta_mixture];
%EM
order = 1:n;
[beta_EM, sigma_EM] = EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, 0, beta_naive, order);
%EMM
order = 1:n;
theta = Choose_theta(n,K);
[beta_EMM, sigma_EMM] =  EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);

error_est(i,j,1) = norm(beta_naive - beta)/b;
error_est(i,j,2) = norm(beta_oracle  - beta)/b;
error_est(i,j,3) = norm(beta_EM  - beta)/b;
error_est(i,j,4) = norm(beta_EMM  - beta)/b;
error_est(i,j,5) = norm(beta_robust  - beta)/b;
error_est(i,j,6) = norm(beta_mixture  - beta)/b;
end
end

save('Hamming_LR.mat');
