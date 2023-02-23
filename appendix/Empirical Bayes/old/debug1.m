addpath('..\functions')
iter = 200; % number of iteration for EM alogrithm
mcmc_steps = 4000; burn_steps = 2000; 
n = 10; d = 1; 
X =  (1:1:10)'; b = 1; sigma = 0.1;
rng(1)
Y = X*b + sigma*randn(n,1);
%plot(X,Y,'*');
Y_P = Y;
Y_P(3) = Y(8);
Y_P(8) = Y(3);
plot(X,Y_P,'*')
P_true = [1 2 8 4 5 6 7 3 9 10];
%-------------------------------
beta_naive = X\Y_P;
%-------------------------------
order = 1:n;
theta = 0;
[beta_EMEB, sigma_EMEB, track1, track2] =  EM_mal_EB(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);

beta_EMEB

b

sort(X)\sort(Y_P)

plot(track1, '*')
plot(track2, '*')

%EM
order = 1:n;
[beta_EM, sigma_EM] = EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, 0, beta_naive, order);

