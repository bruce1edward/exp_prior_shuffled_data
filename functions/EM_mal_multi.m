function [B_hat2, sigma_sq] = EM_mal_multi(Y_permuted, X, iters, mcmc_steps, burn_steps, theta, beta_start, order)
[n,m] = size(Y_permuted);
B_hat2 = beta_start; % initial OLS estimator of B
sigma_sq = norm(Y_permuted - X*B_hat2,'fro')^2/(n*m); % initial OLS estimator of sigma^2
for k = 1 : iters
        Y_hat = X*B_hat2/sigma_sq; % Computing the estimator of Y_hat
        [hat_Pi,order] = mcmc_mex_mal(Y_permuted, Y_hat, order, mcmc_steps, burn_steps, theta);
        hat_Pi = hat_Pi/(mcmc_steps - burn_steps); % take the average of estimating Pi
        B_hat2 = X\(hat_Pi'*Y_permuted);    % Get the Least square solution 
        sigma_sq = norm(Y_permuted - X*B_hat2,'fro')^2/(n*m);
end
end