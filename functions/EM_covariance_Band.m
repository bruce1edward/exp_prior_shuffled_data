function [S_pi] = EM_covariance_Band(x, y, order, mcmc_steps, burn_steps, iters, omega_start, p, q, r)
%{
n = 1000; p = 1; q = 1; sigmax = 1; sigmay = 1; rho = 0.5; k = floor(0.4*n);
[x,y,Pi] = generate_distribution_covariance(n, p, q, sigmax, sigmay, rho, k);
z_naive = [x(Pi,:), y];S_naive = cov(z_naive)*(n-1)/n;
z = [x , y];S_oracle = cov(z)*(n-1)/n;
order = 1:n; mcmc_steps = 8000; burn_steps = 4000;iters = 200;omega_start = S_naive\eye(p+q);
theta = log(n-k);
%}

%{
x = x(Pi,:);iters = iter;omega_start = S_naive\eye(p+q);
%}
n = numel(order);
omega_hat = omega_start; % initial OLS estimator of Omega
XTX = x'*x/n;
YTY = y'*y/n;
for k = 1 : iters
        [hat_Pi,order,count,count_step] = mcmc_covariance_band(x, y, order, mcmc_steps, burn_steps, omega_hat, r);
        yTPx = y'*(hat_Pi*x/count_step)/n;   % compute Y^T Pi X term
 %{     
        yTPx1 = y'*hat_Pi*x/(n*(mcmc_steps - burn_steps))
        xTPty = x'*Pi_y/n;
        cov([Pi_x, y])*(n-1)/n
        Pi_x'*Pi_x/n
 %}
        S_pi = [XTX, yTPx'; yTPx YTY];
        omega_hat = S_pi\eye(p+q);
end
end