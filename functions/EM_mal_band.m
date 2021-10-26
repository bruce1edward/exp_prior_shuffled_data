function [B_hat2, sigma_sq] = EM_mal_band(Y_permuted, X, iters, mcmc_steps, burn_steps, beta_start, order, r)
n = numel(order);%Y_permuted = Y_B; iters = iter;beta_start = beta_naive;
B_hat2 = beta_start; % initial OLS estimator of B
sigma_sq = norm(Y_permuted - X*B_hat2)^2/n; % initial OLS estimator of sigma^2
for k = 1 : iters
        Y_hat = X*B_hat2/sigma_sq; % Computing the estimator of Y_hat
        [hat_Pi, order, count, count_steps] = mcmc_band1(Y_permuted, Y_hat, order, n, mcmc_steps, burn_steps, r);
        hat_Pi = hat_Pi/count_steps;
        Pi_Y = hat_Pi'*Y_permuted;
        B_hat2 = X\Pi_Y;    % Get the Least square solution 
        sigma_sq = norm(Pi_Y - X*B_hat2)^2/n;
end
%norm(B_hat2  - beta)/b
end