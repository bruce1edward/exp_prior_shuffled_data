rng(123)
addpath('..\functions')
iter = 200; % number of iteration for EM alogrithm
mcmc_steps = 4000; burn_steps = 2000; 
%---------------------------------------------------
%n = 1000; d= 5; K = floor(n*0.3); sigma = 1; b = 3;
%[X,Y,Y_P,beta,Pi,inv_Pi] = generate_distribution_sparse(n, d, K, sigma, b);
%---------------------------------------------------
%n = 2000; d= 5;
%r = floor(10); sigma = 1; b = 3;
%[X,Y,Y_P,beta,Pi,average] = generate_distribution_band1(n, d, r, sigma, b);
%----------------------------------------------------
n = 1000; d= 5; r = 50; K = floor(r*0.3); sigma = 1; b = 4; 
[X,Y,Y_P,beta,Pi,blockindex,Q] = generate_distribution_block(n, d, r, K, sigma, b);

%Naive 
beta_naive = X\Y_P;
%Oracle
beta_oracle = X\Y;
%Robust
%beta_robust = robustfit(X_wo,Y_P,'huber');
%EM
order = 1:n;
[beta_EM, sigma_EM] = EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, 0, beta_naive, order);
%EMM
%order = 1:n;
%theta = Choose_theta(n,K);
%[beta_EMM, sigma_EMM] =  EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);
%EMM
%order = 1:n;
%theta = 4.3;
%[beta_EMM1, sigma_EMM1] =  EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);
%EM empirical bayes
order = 1:n;
theta = 0;
[beta_EMEB, sigma_EMEB, track1, track2] =  EM_mal_EB(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);
%end

[norm(beta_naive - beta)/b norm(beta_EM  - beta)/b norm(beta_EMEB - beta)/b ]
