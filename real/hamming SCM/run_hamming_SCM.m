rng('default')
load('../../data/SCM/scm_processed.mat')
k_n = 0.05 : 0.05 : 0.4;
addpath('../../functions')
iter = 200; % number of iteration for EM alogrithm
iters = 100;mcmc_steps = 8000; burn_steps = 4000; nmethod = 3;
error_est = zeros(numel(k_n),iters,nmethod);
%------------------------------------

for i = 1:numel(k_n)
for j = 1:iters
%Sparse Permtuation 
K = floor(n*k_n(i)); 
pi = randperm(K);
Pi = 1:n;
Pi(sort(pi)) = pi;
%Naive 
z_naive = [X(Pi,:) , Y];S_naive = cov(z_naive)*(n-1)/n;R_naive = corrcoef(z_naive);
%Oracle
z_oracle = [X , Y];S_oracle = cov(z_oracle)*(n-1)/n;R_oracle = corrcoef(z_oracle);
%Centered x and y
X_center = X - mean(X,1);
Y_center = Y - mean(Y,1);
%Variance Matrix
var_vector = sqrt(diag(S_oracle)); sigmaxy = zeros(d+m,d+m);
for ii = 1:d+m
    for jj = 1:d+m
        sigmaxy(ii,jj) = var_vector(ii)*var_vector(jj);
    end
end
%EM
order = 1:n;
omega_start = S_naive\eye(d+m); 
S_pi = EM_covariance1(X_center(Pi,:), Y_center, order, mcmc_steps, burn_steps, iter, omega_start, d, m, 0);
%EMM
order = 1:n;
omega_start = S_naive\eye(d+m); 
theta = Choose_theta(n,K);
S_pi_mallows = EM_covariance1(X_center(Pi,:), Y_center, order, mcmc_steps, burn_steps, iter, omega_start, d, m, theta);

error_est(i,j,1) = norm(R_naive - R_oracle, 'fro');
error_est(i,j,2) = norm(S_pi./sigmaxy - R_oracle, 'fro');
error_est(i,j,3) = norm(S_pi_mallows./sigmaxy - R_oracle, 'fro');
end
end