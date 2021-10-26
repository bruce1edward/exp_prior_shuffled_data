function [B_hat2] = EM_mal_Block_GLM1(Y_permuted, X, iters, mcmc_steps, burn_steps, beta_start, g_ID, Y)
%Y_permuted = Y_P; beta_start = beta_naive; X = X1;iters = iter;
uni_g_ID = unique(g_ID);C = cell(numel(uni_g_ID),3);n = numel(Y_permuted);
for g = 1:numel(uni_g_ID)
C{g,1} = find(g_ID == uni_g_ID(g))';C{g,2} = 1:numel(C{g,1});
K_N = sum(Y_permuted(C{g,1}) ~= Y(C{g,1}));
theta = Choose_theta(numel(C{g,1}), K_N);
if isnan(theta)
    C{g,3} = 0;
else
    C{g,3} = theta;
end
end
B_hat2 = beta_start; % initial OLS estimator of B
for k = 1 : iters
        Y_hat = [ones(n,1) X]*B_hat2; % Computing the estimator of Y_hat
        Pi_Y = zeros(n,1);
        for g = 1:numel(uni_g_ID) % g = 20
        [Pi_Y(C{g,1}),C{g,2}] = mcmc_mex_mal_tricks(Y_permuted(C{g,1}), Y_hat(C{g,1}), C{g,2}, mcmc_steps, burn_steps, C{g,3});
        end
        glm = fitglm(X, Pi_Y,'linear','Distribution','poisson');
        B_hat2 = glm.Coefficients.Estimate;    % Get the Least square solution   
end
%norm(beta_EM  - beta_oracle)
%norm(B_hat2 - beta_oracle)
%[norm(beta_naive - beta_oracle), norm(beta_EM  - beta_oracle), norm(B_hat2  - beta_oracle)]
%[dev(beta_naive), dev(beta_EM), dev(B_hat2)]
end