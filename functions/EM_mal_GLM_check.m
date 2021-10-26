function [B_hat2] = EM_mal_GLM_check(Y_permuted, X, iters, mcmc_steps, burn_steps, theta, beta_start, order)
B_hat2 = beta_start; % initial OLS estimator of B
%B_hat2 = beta_naive;Y_permuted = Y_P;%X = [ones(n,1) X];
for k = 1 : iters
        Y_hat = [ones(n,1) X]*B_hat2; % Computing the estimator of Y_hat 
        [hat_Pi,order] = mcmc_mex(Y_permuted, Y_hat, order, mcmc_steps, burn_steps);
        hat_Pi = hat_Pi/(mcmc_steps - burn_steps); % take the average of estimating Pi
        
        %glm = fitglm(X, round(hat_Pi'*Y_permuted),'linear','Distribution','poisson');
        %B_hat3 = glm.Coefficients.Estimate;
        glm1 = fitglm(hat_Pi*X, Y_permuted,'linear','Distribution','poisson');
        B_hat2 = glm1.Coefficients.Estimate;    % Get the Least square solution 
end
end

%sbd = hat_Pi'*Y_permuted;
%avcd = hat_Pi*X;
%X(:,1) = []; \
sbd1  = round(hat_Pi'*Y_permuted);
norm(B_hat2  - beta)/b