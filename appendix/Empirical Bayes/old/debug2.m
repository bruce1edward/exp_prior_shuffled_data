rng(123)
% ORIGINAL: n = 1000; d= 1; K = floor(n*0.3); sigma = 1; b = 1;
n = 100; d= 1; K = floor(n*0.2); sigma = 0.1; b = 1;
[X,Y,Y_P,beta,Pi,inv_Pi] = generate_distribution_sparse(n, d, K, sigma, b);
%Naive 
beta_naive = X\Y_P;

%EM EB
iters = 20;iter = 400; % number of iteration for EM alogrithm
mcmc_steps = 8000; burn_steps = 4000; 
log_psi = @(gamma) exp(gamma) - 1 - gamma*n;
track_theta = zeros(iters,1);
track_sum = zeros(iters,1);
Y_permuted = Y_P; iters = iter; beta_start = beta_naive; order = 1:n; theta = log(n*.15);%Choose_theta(n,K);
n = numel(order);
B_hat2 = beta_start; % initial OLS estimator of B
sigma_sq = norm(Y_permuted - X*B_hat2)^2/n; % initial OLS estimator of sigma^2
for k = 1 : iters
        Y_hat = X*B_hat2/sigma_sq; % Computing the estimator of Y_hat
        [hat_Pi,order] = mcmc_mex_mal(Y_permuted, Y_hat, order, mcmc_steps, burn_steps, theta);
        hat_Pi = hat_Pi/(mcmc_steps - burn_steps);
        Pi_Y = hat_Pi*Y_permuted;
        B_hat2 = X\(hat_Pi*Y_permuted);    % Get the Least square solution 
        sigma_sq = mean(Y_P.^2) + mean((X*B_hat2).^2) - 2 * mean(Pi_Y .* (X*B_hat2));%/% norm(Pi_Y - X*B_hat2)^2/n;
        %fun = @(gamma) n*log_psi(gamma) + gamma*sum(diag(hat_Pi));
        theta = log(sum(diag(hat_Pi)));%fminbnd(fun,0,2*log(n));
        track_theta(k) = theta;
        track_sum(k) = sum(diag(hat_Pi));
end

gr = 0:0.01:(2*log(n));
plot(gr, fun(gr))

%Optimal Gamma in terms of hamming distance
%Sorting
sort_pi = E_Pi1(Y_P,X*beta);
%sum(sort(Y_P).*sort(X))
%sum(Y_P(sort_pi).*X)
%MAP
theta = 0:0.1:3*log(n);
track_map = zeros(numel(theta),1);
for ii = 1:numel(theta)
%ii=1;
[assignments,P] = MAP_est_multi(n, Y_P, X, theta(ii), sigma);
track_map(ii) = norm(Y - Y_P(assignments))/norm(Y);
end


figure 
hold on
p1 = plot(theta, track_map,'ro');
p3 = yline(norm(Y - Y_P(sort_pi))/norm(Y),'LineWidth',1.5);
p4 = xline(log(n),'LineWidth',1.5);
p5 = xline(Choose_theta(n,K),'LineWidth',1.5);
%title(['n = ' num2str(n) ' K = ' num2str(K/n) ' SNR= ' num2str(b^2/sigma_sq)])
xlabel('\theta')
ylabel('$||\hat{\Pi}^{MAP} Y - \Pi^{*} Y||_{2}/||Y||_{2}$','interpreter','latex')
legend([p1 p3 p4 p5], {'MAP', 'sorting','log(n) - Empirical Bayes', 'optimal via concentration'}, 'location', 'best');
hold off

export_fig('eb5.pdf')

[a,b] = min(track_map);
theta(b)/log(n)
