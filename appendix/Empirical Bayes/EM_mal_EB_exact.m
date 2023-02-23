function [B_hat2, sigma_sq, track_theta, track_sum, track_fval, track_exitflag] = EM_mal_EB_exact(Y_permuted, X, iters, n, theta, b)
track_theta = zeros(iters,1);
track_sum = zeros(iters,1);
track_fval = zeros(iters,1);
track_exitflag = zeros(iters,1);
B_hat2 = b; % initial OLS estimator of B
sigma_sq = norm(Y_permuted - X*B_hat2)^2/n; % initial OLS estimator of sigma^2
for k = 1 : iters
        %[hat_Pi,order] = mcmc_mex_mal(Y_permuted, Y_hat, order, mcmc_steps, burn_steps, theta);
        %hat_Pi = hat_Pi/(mcmc_steps - burn_steps);
        hat_Pi = exact_EM(Y_permuted, X, B_hat2, sigma_sq, n, theta);
        fun = @(gamma) exp(gamma) - 1 - gamma*sum(diag(hat_Pi));
        [theta,fval,exitflag] = fminbnd(fun,0,2*log(n));
        track_fval(k) = fval;
        track_exitflag(k) = exitflag;
        track_theta(k) = theta;
        track_sum(k) = sum(diag(hat_Pi));
        Pi_Y = hat_Pi*Y_permuted;
        B_hat2 = X\Pi_Y;    % Get the Least square solution 
        sigma_sq = mean(Y_permuted.^2) + mean((X*B_hat2).^2) - 2 * mean(Pi_Y .* (X*B_hat2));
end
end