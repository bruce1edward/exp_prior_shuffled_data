rng('default')
addpath('..\..\functions')
K_n = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4];
iters = 100;
error_est = zeros(iters,numel(K_n),7);
r_square = zeros(iters,numel(K_n),7);
%---------------------------------Data(Sparse Permutating)
data = readtable('../../data/CPS/fakedata.csv','ReadVariableNames',false);
Y = data(2:end,2);X = data(2:end,3:8);
Y = Y{:,:};X = X{:,:};
Y = cellfun(@str2num,Y);X = cellfun(@str2num,X);
%hot encoding
cat = unique(X(:,5));
[n,d] = size(X);
for i = 1:numel(cat)-1
    X(:,d+i) = X(:,5) == cat(i+1);
end
temp1 = X(:,6);
X(:,5:6) = [];
X = [X temp1];
X = [ones(n,1) X];
X_wo = X(:,2:end);
X_centered = X_wo - repmat(mean(X_wo,1), [n 1]);
%---------------------------------------------------------------
for i = 1:iters
for j = 1:numel(K_n)
pi = randperm(floor(K_n(j)*numel(Y)));Pi = 1:numel(Y);Pi(sort(pi)) = pi;
Y_permuted = Y(Pi);
[n,d] = size(X);
[n,m] = size(Y);
%Naive
beta_naive = X\Y_permuted;
sigma_sq_naive = norm(Y_permuted - X*beta_naive)^2/n;
r_sq_naive = 1 - norm(Y - X*beta_naive)^2/norm(Y - mean(Y))^2;
residual_variance_naive = norm(Y - X*beta_naive)^2/(n-3);
%Oracle
beta_oracle = X\Y;
sigma_sq_oracle = norm(Y - X*beta_oracle)^2/n;
r_sq_oracle = 1 - norm(Y - X*beta_oracle)^2/norm(Y - mean(Y))^2;
residual_variance_oracle = norm(Y - X*beta_oracle)^2/(n-3);
%Robust
beta_robust = robustfit(X_wo,Y_permuted,'huber');
%EM
iter = 200; mcmc_steps = 4000; burn_steps = 2000; order = 1:n; 
[beta_EM, sigma_sq_EM] = EM_mal_tricks(Y_permuted, X, iter, mcmc_steps, burn_steps, 0, beta_naive, order);
%EMM
iter = 200; mcmc_steps = 4000; burn_steps = 2000; order = 1:n; theta = Choose_theta(n,floor(K_n(j)*n));
[beta_EMM, sigma_sq_EMM] = EM_mal_tricks(Y_permuted, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);
%Mixture
Y_permuted1 = Y_permuted - mean(Y_permuted); control = "robust";
[beta_mixture, alphahat, sigmahat] = fit_mixture(X_centered, Y_permuted1, control);
beta_mixture = [mean(Y_permuted - X_wo*beta_mixture) ; beta_mixture];
%LL
Q = (1 - K_n(j) - K_n(j)/(n-1)) * eye(n) + (K_n(j)/(n-1)) * ones(n);
QX = Q*X;
beta_LL = QX\Y_permuted;
%Chamber
beta_C = (X'*Q*X)\(X'*Y_permuted);

error_est(i,j,1) = norm(beta_naive - beta_oracle)/norm(beta_oracle);
error_est(i,j,2) = norm(beta_robust  - beta_oracle)/norm(beta_oracle);
error_est(i,j,3) = norm(beta_mixture  - beta_oracle)/norm(beta_oracle);
error_est(i,j,4) = norm(beta_EM  - beta_oracle)/norm(beta_oracle);
error_est(i,j,5) = norm(beta_EMM  - beta_oracle)/norm(beta_oracle);
error_est(i,j,6) = norm(beta_LL  - beta_oracle)/norm(beta_oracle);
error_est(i,j,7) = norm(beta_C  - beta_oracle)/norm(beta_oracle);
end
end

