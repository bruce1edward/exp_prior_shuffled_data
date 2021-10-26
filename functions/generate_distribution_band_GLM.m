function [X,Y_B,beta,Pi,blockindex] = generate_distribution_band_GLM(n, d, r, b)
%rng(3) % for reproducibility
X = randn(n,d);
X_w_intercept = [ones(n,1) X];
%X = rand(n,d-1)*2*sqrt(3) - sqrt(3);
beta = normrnd(0,1,[d + 1,1]);
beta = beta/norm(beta);
beta = beta*b;
%Band Permuted Data
Pi = zeros(n,1);
blockindex = zeros(n,1);
for i = 1:n/r
   Pi(1 + (i-1)*r: i*r) = randperm(r) + (i-1)*r;
   blockindex(1+ (i-1)*r: i*r) = i;
end
mu = exp(X_w_intercept(Pi,:)*beta);
Y_B  = poissrnd(mu);
end