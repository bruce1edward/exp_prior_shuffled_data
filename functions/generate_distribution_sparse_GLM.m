function [X,Y_P,beta,Pi,inv_Pi] = generate_distribution_sparse_GLM(n, d, K, b)
%rng(1) % for reproducibility
X = randn(n,d);
X_w_intercept = [ones(n,1) X];
%X = rand(n,d-1)*2*sqrt(3) - sqrt(3);
beta = normrnd(0,1,[d + 1,1]);
beta = beta/norm(beta);
beta = beta*b;
%Sparse Permuted Data
pi = randperm(K);
Pi = 1:n;
Pi(sort(pi)) = pi;
%sum(Pi~=1:n)
mu = exp(X_w_intercept(Pi,:)*beta);
Y_P = poissrnd(mu);
[abc, inv_Pi] = sort(Pi);
end

