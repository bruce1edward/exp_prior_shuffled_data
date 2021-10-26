function [B_hat2, sigma_sq] = EM_mal_Block(Y_permuted, X, iters, mcmc_steps, burn_steps, theta, beta_start, g_ID)
%beta_start = beta_naive;order = 1:n;theta = 0;mcmc_steps = 2000;
%burn_steps = 500;iters = 200;g_ID = blockindex;Y_permuted = Y_P;
uni_g_ID = unique(g_ID);C = cell(numel(uni_g_ID),2);n = numel(Y_permuted);
for g = 1:numel(uni_g_ID)
C{g,1} = find(g_ID == uni_g_ID(g))';C{g,2} = 1:numel(C{g,1});
end
B_hat2 = beta_start; % initial OLS estimator of B
sigma_sq = norm(Y_permuted - X*B_hat2)^2/n; % initial OLS estimator of sigma^2
for k = 1 : iters
        Y_hat = X*B_hat2/sigma_sq; % Computing the estimator of Y_hat
        Pi_Y = zeros(n,1);
        for g = 1:numel(uni_g_ID)
        [Pi_Y(C{g,1}),C{g,2}] = mcmc_mex_mal_tricks(Y_permuted(C{g,1}), Y_hat(C{g,1}), C{g,2}, mcmc_steps, burn_steps, theta);
        end
        B_hat2 = X\Pi_Y;    % Get the Least square solution 
        sigma_sq = norm(Pi_Y - X*B_hat2)^2/n;
end
end