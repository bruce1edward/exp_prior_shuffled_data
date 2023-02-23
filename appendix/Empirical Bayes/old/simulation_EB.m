addpath('..\functions')
iters = 20;iter = 400; % number of iteration for EM alogrithm
mcmc_steps = 8000; burn_steps = 4000; 
n = 1000; d= 20; K = floor(n*0.3); sigma = 1; b = 3;
[X,Y,Y_P,beta,Pi,inv_Pi] = generate_distribution_sparse(n, d, K, sigma, b);
X_wo = X(:,2:end);X_centered = X_wo - repmat(mean(X_wo,1), [n 1]);
%Naive 
beta_naive = X\Y_P;
%EMH
track_theta = [linspace(0,log(n-K),10), linspace(log(n-K),log(n),100)];
track_mse = zeros(numel(track_theta), 1);
for i = 1:numel(track_theta)
order = 1:n;
theta = track_theta(i);
[beta_EMM, sigma_EMM] =  EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);
track_mse(i) = norm(beta_EMM  - beta)/b;
end

plot(track_theta(1:10), track_mse(1:10), '*')
xlabel('\gamma', 'FontSize', 14)
ylabel('mse(\gamma)', 'FontSize', 14)
export_fig('eb3.pdf')

plot(track_theta(11:110), track_mse(11:110), '*')
xlabel('\gamma', 'FontSize', 14)
ylabel('mse(\gamma)', 'FontSize', 14)
p1 = xline(log(n-K));
p2 = xline( Choose_theta(n,K),'--r');
p3 = xline(log(n), '--b');
lh = legend([p1 p2 p3], {'log(n-K)','optimal \gamma','log(n)'}, 'location', [0.3416 0.6282 0.3661 0.1833]);
export_fig('eb4.pdf')