function [B_hat2] = EM_mal_GLM(n, Y_permuted, X, iters, mcmc_steps, burn_steps, theta, beta_start, order)
B_hat2 = beta_start; % initial OLS estimator of B
%B_hat2 = beta_naive;Y_permuted = Y_P;theta = 0;
for k = 1 : iters
        Y_hat = [ones(n,1) X]*B_hat2; % Computing the estimator of Y_hat
        [hat_Pi,order] = mcmc_mex_mal(Y_permuted, Y_hat, order, mcmc_steps, burn_steps, theta);
        hat_Pi = hat_Pi/(mcmc_steps - burn_steps); % take the average of estimating Pi
        glm = fitglm(X, hat_Pi'*Y_permuted,'linear','Distribution','poisson');
        B_hat2 = glm.Coefficients.Estimate;    % Get the Least square solution 
end
end