rng('default')
addpath('..\..\functions')
k_n = [0.2 0.25 0.3 0.35 0.4 0.45 0.5];nmethod = 5;
p = 5; q = 5;sigmax = 1; sigmay = 1; rho1 = 0.8; iter = 400; % number of iteration for EM alogrithm
iters = 100;
mcmc_steps = 8000; burn_steps = 4000; 
error_est = zeros(numel(k_n),iters,nmethod);
r_square = zeros(numel(k_n),iters,nmethod);

for i = 1:numel(k_n)
for j = 1:iters
n = 1000; K = floor(n*k_n(i)); 
[x,y,Pi,sigma] = generate_distribution_covariance(n, p, q, sigmax, sigmay, rho1, K);
%Naive 
z_naive = [x(Pi,:) , y];S_naive = cov(z_naive)*(n-1)/n;
%Oracle
z_oracle = [x , y];S_oracle = cov(z_oracle)*(n-1)/n;
%EM
order = 1:n;
omega_start = S_naive\eye(p+q); 
S_pi = EM_covariance1(x(Pi,:), y, order, mcmc_steps, burn_steps, iter, omega_start, p, q, 0);
%EMM
order = 1:n;
theta = Choose_theta(n,K)/sqrt(p+q);
omega_start = S_naive\eye(p+q); 
S_pi_mallows = EM_covariance1(x(Pi,:), y, order, mcmc_steps, burn_steps, iter, omega_start, p, q, theta);
%robust
sig_robust = robustcov(z_naive);

error_est(i,j,1) = norm(S_naive - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,2) = norm(S_oracle - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,3) = norm(S_pi - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,4) = norm(S_pi_mallows - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,5) = norm(sig_robust - sigma, 'fro')/norm(sigma, 'fro');
end
end

save('Hamming_MVN.mat');