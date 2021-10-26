function [B_hat2, sigma_sq] = EM_mal_tricks(Y_permuted, X, iters, mcmc_steps, burn_steps, theta, beta_start, order)
n = numel(order);
B_hat2 = beta_start; % initial OLS estimator of B
sigma_sq = norm(Y_permuted - X*B_hat2)^2/n; % initial OLS estimator of sigma^2
for k = 1 : iters
        Y_hat = X*B_hat2/sigma_sq; % Computing the estimator of Y_hat
        [Pi_Y,order] = mcmc_mex_mal_tricks(Y_permuted, Y_hat, order, mcmc_steps, burn_steps, theta);
        B_hat2 = X\Pi_Y;    % Get the Least square solution 
        sigma_sq = norm(Pi_Y - X*B_hat2)^2/n;
end
end