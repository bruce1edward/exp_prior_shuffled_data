function [X,Y,Y_B,beta,Pi,blockindex] = generate_distribution_band(n, d, r, sigma, b)
%rng(3) % for reproducibility
X = randn(n,d-1);
X = [ones(n,1) X];
beta = normrnd(0,1,[d,1]);
beta = beta/norm(beta);
beta = beta*b;
%Band Permuted Data
Pi = zeros(n,1);
blockindex = zeros(n,1);
for i = 1:n/r
   Pi(1 + (i-1)*r: i*r) = randperm(r) + (i-1)*r;
   blockindex(1+ (i-1)*r: i*r) = i;
end
[abc, inv_Pi] = sort(Pi);
Y = X*beta + normrnd(0,sigma,[n,1]);
Y_B = Y(Pi);
end