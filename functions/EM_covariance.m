function [S_pi] = EM_covariance(x, y, order, mcmc_steps, burn_steps, iters, omega_start, p, q, theta)
%{
n = 1000; p = 1; q = 1; sigmax = 1; sigmay = 1; rho = 0.5; k = floor(0.4*n);
[x,y,Pi] = generate_distribution_covariance(n, p, q, sigmax, sigmay, rho, k);
z_naive = [x(Pi,:), y];S_naive = cov(z_naive)*(n-1)/n;
z = [x , y];S_oracle = cov(z)*(n-1)/n;
order = 1:n; mcmc_steps = 8000; burn_steps = 4000;iters = 200;omega_start = S_naive\eye(p+q);
theta = log(n-k);
%}

%{
%sum(Pi ~= 1:n)/n
x = x_center(Pi,:);
y = y_center;
iters = iter;
omega_start = S_naive\eye(p+q);
order = 1:n;
%theta = 0;
theta = Choose_theta(n,K);
%
%}
n = numel(order);
omega_hat = omega_start; % initial OLS estimator of Omega
XTX = x'*x/n;
YTY = y'*y/n;

for k = 1 : iters
        [Pi_x,order,count] = mcmc_covariance(x, y, order, mcmc_steps, burn_steps, omega_hat, theta);
        yTPx = y'*Pi_x/n;   % compute Y^T Pi X term
 %{     
        yTPx1 = y'*hat_Pi*x/(n*(mcmc_steps - burn_steps))
        xTPty = x'*Pi_y/n;
        cov([Pi_x, y])*(n-1)/n
        Pi_x'*Pi_x/n
 %}
        S_pi = [XTX, yTPx'; yTPx YTY];
        omega_hat = S_pi\eye(p+q);
end
%norm(S_pi - S_oracle, 'fro')
end