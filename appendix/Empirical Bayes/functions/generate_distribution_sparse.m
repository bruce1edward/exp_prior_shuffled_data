function [X,Y,Y_P,beta,Pi,inv_Pi] = generate_distribution_sparse(n, d, K, sigma, b)
%rng(1) % for reproducibility
%X = randn(n,d-1);
%X = [ones(n,1) X];
X = randn(n,d);
beta = normrnd(0,1,[d,1]);
beta = beta/norm(beta);
beta = beta*b;
%Sparse Permuted Data
pi = randperm(K);
Pi = 1:n;
Pi(sort(pi)) = pi;
%sum(Pi~=1:n)
Y = X*beta + normrnd(0,sigma,[n,1]);
Y_P = Y(Pi);
[abc, inv_Pi] = sort(Pi);
end