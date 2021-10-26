rng('default')
addpath('..\..\functions')
k_n = [0.2 0.25 0.3 0.35 0.4 0.45 0.5];nmethod = 5;
p = 5; q = 5;
iters = 100;iter = 400; % number of iteration for EM alogrithm
mcmc_steps = 4000; burn_steps = 2000; 
error_est = zeros(numel(k_n),iters,nmethod);
r_square = zeros(numel(k_n),iters,nmethod);

for i = 1:numel(k_n)
for j = 1:iters
n = 1000; r = 50; K = floor(r*k_n(i)); sigmax = 1; sigmay = 1; rho1 = 0.8; 
[x,y,Pi,sigma,blockindex,Q] = generate_distribution_covariance_block(n, p, q, sigmax, sigmay, rho1, K, r);
%Naive 
z_naive = [x(Pi,:) , y];S_naive = cov(z_naive)*(n-1)/n;
%Oracle
z_oracle = [x , y];S_oracle = cov(z_oracle)*(n-1)/n;
%LL-C
QX = Q*x; XTX = x'*x/n; YTY = y'*y/n; yTPx = y'*(QX - mean(QX))/n;
S_LL = [XTX, yTPx'; yTPx YTY];
%C
%EM
omega_start = S_naive\eye(p+q); 
S_pi = EM_covariance_Block(n, x(Pi,:), y, mcmc_steps, burn_steps, iter, omega_start, p, q, 0, blockindex);
%EMM
theta = Choose_theta(r, K);
omega_start = S_naive\eye(p+q); 
S_pi_mallows = EM_covariance_Block(n, x(Pi,:), y, mcmc_steps, burn_steps, iter, omega_start, p, q, theta, blockindex);

error_est(i,j,1) = norm(S_naive - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,2) = norm(S_oracle - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,3) = norm(S_pi - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,4) = norm(S_pi_mallows - sigma, 'fro')/norm(sigma, 'fro');
error_est(i,j,5) = norm(S_LL - sigma, 'fro')/norm(sigma, 'fro');
end
end

save('Block_MVN.mat');