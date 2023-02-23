rng('default')
addpath('functions') 
addpath('data')
load('Italian_survey_data_linkage.mat');X_link = X;X = X(:,1);Y_link = Y;Y = Y(:,1);
n = numel(X);X = [ones(n,1) X];K_n = 0.4;
track_ree = zeros(100, 4);
%sum(Pi~=1:n)
for i = 1:100
pi_ = randperm(floor(K_n*numel(Y)));Pi1 = 1:numel(Y);Pi1(sort(pi_)) = pi_;
Y_permuted = Y(Pi1);
%Naive
beta_naive = X\Y_permuted;
sigma_sq_naive = norm(Y_permuted - X*beta_naive)^2/n;
%Oracle
beta_oracle = X\Y;
sigma_sq_oracle = norm(Y - X*beta_oracle)^2/n;
%----------------------------------------------
iter = 200;mcmc_steps = 4000;burn_steps = 2000;
%EM
order = 1:n;
[beta_EM, sigma_EM] = EM_mal_tricks(Y_permuted, X, iter, mcmc_steps, burn_steps, 0, beta_naive, order);
%EMM
order = 1:n;
theta = Choose_theta(n,floor(K_n*n));
[beta_EMM, sigma_EMM] =  EM_mal_tricks(Y_permuted, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);
%EM Empirical Bayes
order = 1:n; theta = Choose_theta(n,floor(K_n*n));%log(n);
[beta_EMMB, sigma_EMMB, track_theta] =  EM_mal_EB(Y_permuted, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);

track_ree(i,1) = norm(beta_naive - beta_oracle);
track_ree(i,2) =norm(beta_EM - beta_oracle);
track_ree(i,3) =norm(beta_EMM - beta_oracle);
track_ree(i,4) =norm(beta_EMMB - beta_oracle);
end
