function [S_pi] = EM_covariance_Block(n, x, y, mcmc_steps, burn_steps, iters, omega_start, p, q, theta, g_ID)
uni_g_ID = unique(g_ID);C = cell(numel(uni_g_ID),2);
for g = 1:numel(uni_g_ID)
C{g,1} = find(g_ID == uni_g_ID(g))';C{g,2} = 1:numel(C{g,1});
end
omega_hat = omega_start; % initial OLS estimator of Omega
XTX = x'*x/n;
YTY = y'*y/n;
for k = 1 : iters
        Pi_x = zeros(n,p);
        for g = 1:numel(uni_g_ID)
        [Pi_x(C{g,1},:),C{g,2},count] = mcmc_covariance(x(C{g,1},:), y(C{g,1},:), C{g,2}, mcmc_steps, burn_steps, omega_hat, theta);
        end
        yTPx = y'*Pi_x/n;   % compute Y^T Pi X term
        S_pi = [XTX, yTPx'; yTPx YTY];
        omega_hat = S_pi\eye(p+q);
end
end