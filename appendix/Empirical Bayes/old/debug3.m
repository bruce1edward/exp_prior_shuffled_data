rng(123)
n = 10; d= 1; K = floor(n*0.3); sigma = 1; b = 1;
[X,Y,Y_P,beta,Pi,inv_Pi] = generate_distribution_sparse(n, d, K, sigma, b);
%Naive 
beta_naive = X\Y_P;
%Exact Expectation
pi = perms(1:n);
order = 1:n;
spi = pi(sum(pi ~= order,2) <= K,:);
[nsp,d1] = size(spi);
Pi_id = eye(n);
p_density = zeros(nsp,1);
hat_SPi = zeros(n,n);
for j = 1 : nsp 
p_density(j) = dot(Y_P,X(spi(j,:),:)*beta)/sigma;
end
p_density =  p_density - max(p_density);
p_density = softmax(p_density);

for l = 1 : nsp
hat_SPi = hat_SPi + Pi_id(spi(l,:),:)*p_density(l);
end

%Empirical Bayes with exact expectation
iters = 20;iter = 400; % number of iteration for EM alogrithm
mcmc_steps = 8000; burn_steps = 4000; 
log_psi = @(gamma) exp(gamma) - 1 - gamma*n;
track_theta = zeros(iters,1);
track_sum = zeros(iters,1);
Y_permuted = Y_P; iters = iter; beta_start = beta_naive; order = 1:n; theta = log(n);
n = numel(order);
B_hat2 = beta_start; % initial OLS estimator of B
sigma_sq = norm(Y_permuted - X*B_hat2)^2/n; % initial OLS estimator of sigma^2
for k = 1 : iters
        Y_hat = X*B_hat2/sigma_sq; % Computing the estimator of Y_hat
        %[hat_Pi,order] = mcmc_mex_mal(Y_permuted, Y_hat, order, mcmc_steps, burn_steps, theta);
        %hat_Pi = hat_Pi/(mcmc_steps - burn_steps);
        p_density = zeros(nsp,1);
        hat_SPi = zeros(n,n);
        for j = 1 : nsp 
        p_density(j) = dot(Y_hat,X(spi(j,:),:)*B_hat2)/sigma_sq - theta*(sum(spi(j,:)~=1:n));
        end
        p_density =  p_density - max(p_density);
        p_density = softmax(p_density);

        for l = 1 : nsp
        hat_Pi = hat_SPi + Pi_id(spi(l,:),:)*p_density(l);
        end
        
        Pi_Y = hat_Pi*Y_permuted;
        B_hat2 = X\(hat_Pi*Y_permuted);    % Get the Least square solution 
        sigma_sq = norm(Pi_Y - X*B_hat2)^2/n;
        fun = @(gamma) n*log_psi(gamma) + gamma*sum(diag(hat_Pi));
        theta = fminbnd(fun,0,2*log(n));
        track_theta(k) = theta;
        track_sum(k) = sum(diag(hat_Pi));
end
log(n)
norm(B_hat2 - beta)

%EMH
order = 1:n;
theta = Choose_theta(n,K);
[beta_EMM, sigma_EMM] =  EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);
norm(beta_EMM - beta)