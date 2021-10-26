rng('default')
addpath('..\..\functions')
r_n = 3:1:10;nmethod = 5;
p = 5; q = 5;sigmax = 1; sigmay = 1; rho1 = 0.8; iter = 400; % number of iteration for EM alogrithm
iters = 100;
mcmc_steps = 8000; burn_steps = 4000; 
error_est = zeros(numel(r_n),iters,nmethod);
for i = 1:numel(r_n)
for j = 1:iters
n = 1000; r = floor(r_n(i)); 
[x,y,Pi,sigma,average] = generate_distribution_covariance_band1(n, p, q, sigmax, sigmay, rho1, r);
%Naive 
z_naive = [x(Pi,:) , y];S_naive = cov(z_naive)*(n-1)/n;
%Oracle
z_oracle = [x , y];S_oracle = cov(z_oracle)*(n-1)/n;
%EM
order = 1:n;
omega_start = S_naive\eye(p+q); 
S_pi = EM_covariance(x(Pi,:), y, order, mcmc_steps, burn_steps, iter, omega_start, p, q, 0);
%EMM
order = 1:n;
omega_start = S_naive\eye(p+q); 
S_pi_mallows = EM_covariance_Band(x(Pi,:), y, order, mcmc_steps, burn_steps, iter, omega_start, p, q, r);
%average
[n__,d] = size(average);
S_average = cov(average)*(n__-1)/n__;

error_est(i,j,1) = norm(S_naive - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,2) = norm(S_oracle - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,3) = norm(S_pi - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,4) = norm(S_pi_mallows - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,5) = norm(S_average - sigma, 'fro')/norm(sigma, 'fro');
end
end
save('Local_MVN.mat');