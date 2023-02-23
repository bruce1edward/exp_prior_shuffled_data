addpath('..\functions')
n = 10; d = 1; 
X = (1:1:n)'; b = 1; sigma = 0.1;
rng(1)
Y = X*b + sigma*randn(n,1);
%plot(X,Y,'*');
P_true = 1:n;
P = 1:n;
for i = 1:2
  P_true(i) = P(n - i);
  P_true(n-i) = P(i);
end
%sum(P_true ~= P)
Y_P = Y(P_true);
%plot(X,Y_P,'*')
%----------------
%Gap induced by the transition from full integration to Monte Carlo Integration (for n = 10, first)
beta_naive = X\Y_P;
sigma_sq = norm(Y_P - X*beta_naive)^2/n;
%Full integration 
hat_Pi_F = exact_EM(Y_P, X, beta_naive, sigma_sq, n, log(n));
%Monte Carlo integration
order = 1:n; theta = log(n); mcmc_steps = 8000; burn_steps = 4000;
Y_hat = X*beta_naive/sigma_sq; % Computing the estimator of Y_hat
[hat_Pi,order] = mcmc_mex_mal(Y_P, Y_hat, order, mcmc_steps, burn_steps, theta);
hat_Pi = hat_Pi/(mcmc_steps - burn_steps);

norm(hat_Pi_F - hat_Pi, 'F')

iter = 2;order = 1:n; theta = log(n);
%EM Exact
[beta_EMMBE, sigma_EMMBE, track_thetaE, track_sumE, track_fval, track_exitflag] = EM_mal_EB_exact(Y_P, X, iter, n, theta, beta_naive);
%EM MH
[beta_EMMB, sigma_EMMB, track_theta, track_sum] =  EM_mal_EB(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);

[beta_EMMBE beta_EMMB]
[sigma_EMMBE sigma_EMMB]
plot(track_thetaE, track_theta, '*')

%Extending the n = 10 example: Increasing n = 20, n = 40, n = 80 and #mismatches, and see what happens as n is increased


%Are there are any issues in the optimization of gamma via the fminbnd function? How is the fmindbnd function initialized (starting point)
fun = @(gamma) exp(gamma) - 1 - gamma*sum(diag(hat_Pi_F));
[theta,fval,exitflag] = fminbnd(fun,0,2*log(n));
theta
fval
exitflag
output