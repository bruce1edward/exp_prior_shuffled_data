rng(123)
addpath('functions')
iter = 200; % number of iteration for EM alogrithm
mcmc_steps = 8000; burn_steps = 4000; 
%---------------------------------------------------
n = 200; d= 1; K = floor(n*0.3); sigma = 1; b = 3;
[X,Y,Y_P,beta,Pi,inv_Pi] = generate_distribution_sparse(n, d, K, sigma, b);
X_wo = X(:,2:end);X_centered = X_wo - repmat(mean(X_wo,1), [n 1]);
%---------------------------------------------------
%Naive 
beta_naive = X\Y_P;
%Oracle
beta_oracle = X\Y;
%EM
order = 1:n;
[beta_EM, sigma_EM] = EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, 0, beta_naive, order);
%EM empirical bayes
order = 1:n;
theta = log(n);
[beta_EMEB, sigma_EMEB, track1, track2] =  EM_mal_EB(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);

hold on
plot(track1,'*');
p1 = yline(log(n), 'r');
p2 = yline(Choose_theta(n,K), 'k');
xlabel('EM iteration', 'FontSize', 14)
ylabel('Empirical Bayes \gamma', 'FontSize', 14)
lh = legend([p1 p2], {'log(n)','concentration \gamma'}, 'location', 'best');
lh.NumColumns = 1;
lh.FontSize = 10;
lh.ItemTokenSize = [30 70];  
legend('boxoff')
export_fig('eb1.pdf')

hold on
plot(track2,'*')
xlabel('EM iteration', 'FontSize', 14)
ylabel('sum(diag(\Pi))', 'FontSize', 14)
%export_fig('eb2.pdf')

hold on
plot(log(track2) - track1,'*')
xlabel('EM iteration', 'FontSize', 14)
ylabel('log(sum(diag(\Pi))) - gamma', 'FontSize', 14)


[norm(beta_naive - beta)/b norm(beta_EM  - beta)/b norm(beta_EMEB - beta)/b norm(beta_EMM  - beta)/b]
%norm(beta_EMM1  - beta)/b

theta = 0:0.1:3*log(n);
track_map = zeros(numel(theta),1);
for ii = 1:numel(theta)
order = 1:n;
[beta_EMM1, sigma_EMM1] =  EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, theta(ii), beta_naive, order);
track_map(ii) = norm(beta_EMM1  - beta)/b;
end

figure 
hold on
p1 = plot(theta, track_map,'ro');
p3 = xline(mean(track1), 'g');
p4 = xline(log(n), 'k');
p5 = xline(Choose_theta(n,K), 'b');
%title(['n = ' num2str(n) ' K = ' num2str(K/n) ' SNR= ' num2str(b^2/sigma_sq)])
xlabel('\gamma')
ylabel('Relative Estimation Error','interpreter','latex')
legend([p3 p4 p5], {'empirical bayes \gamma', 'log(n)', '\gamma via concentration inequality'}, 'location', 'northeast');
hold off

export_fig('eb2.pdf')
