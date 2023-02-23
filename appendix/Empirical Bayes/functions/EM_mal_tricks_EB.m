function [B_hat2, sigma_sq] = EM_mal_tricks_EB(Y_permuted, X, iters, mcmc_steps, burn_steps, theta, beta_start, order)
%Y_permuted = Y_P; iters = iter; beta_start = beta_naive; order = 1:n; theta = Choose_theta(n,K);
n = numel(order);
B_hat2 = beta_start; % initial OLS estimator of B
sigma_sq = norm(Y_permuted - X*B_hat2)^2/n; % initial OLS estimator of sigma^2
log_psi = @(gamma) exp(gamma) - 1 - gamma*n;
for k = 1 : iters
        Y_hat = X*B_hat2/sigma_sq; % Computing the estimator of Y_hat
        [hat_Pi,order] = mcmc_mex_mal(Y_permuted, Y_hat, order, mcmc_steps, burn_steps, theta);
        hat_Pi = hat_Pi/(mcmc_steps - burn_steps);
        Pi_Y = hat_Pi*Y_permuted;
        B_hat2 = X\(hat_Pi*Y_permuted);    % Get the Least square solution 
        sigma_sq = mean(Y_permuted.^2) + mean((X*B_hat2).^2) - 2 * mean(Pi_Y .* (X*B_hat2));%norm(Pi_Y - X*B_hat2)^2/n;
        %fun = @(gamma) n*log_psi(gamma) + gamma*sum(diag(hat_Pi));
        theta = log(sum(diag(hat_Pi)));%fminbnd(fun,0,2*log(n));
        %track_theta(k) = theta;
        %track_sum(k) = sum(diag(hat_Pi));
end