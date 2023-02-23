function [B_hat2, sigma_sq, track_theta, track_sum] = EM_mal_EB(Y_permuted, X, iters, mcmc_steps, burn_steps, theta, beta_start, order)
n = numel(order);
track_theta = zeros(iters,1);
track_sum = zeros(iters,1);
B_hat2 = beta_start; % initial OLS estimator of B
sigma_sq = norm(Y_permuted - X*B_hat2)^2/n; % initial OLS estimator of sigma^2
for k = 1 : iters
        Y_hat = X*B_hat2/sigma_sq; % Computing the estimator of Y_hat
        [hat_Pi,order] = mcmc_mex_mal(Y_permuted, Y_hat, order, mcmc_steps, burn_steps, theta);
        hat_Pi = hat_Pi/(mcmc_steps - burn_steps);
        %fun = @(gamma) exp(gamma) - 1 - gamma*sum(diag(hat_Pi));
        theta = log(sum(diag(hat_Pi)));%%fminbnd(fun,0,2*log(n));
        track_theta(k) = theta;
        track_sum(k) = sum(diag(hat_Pi));%log(sum(diag(hat_Pi)));
        Pi_Y = hat_Pi'*Y_permuted;
        B_hat2 = X\Pi_Y;    % Get the Least square solution 
        sigma_sq = mean(Y_permuted.^2) + mean((X*B_hat2).^2) - 2 * mean(Pi_Y .* (X*B_hat2));%norm(Pi_Y - X*B_hat2)^2/n;
end
end