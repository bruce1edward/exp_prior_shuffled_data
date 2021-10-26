function [X,Y,Y_B,beta,Pi,average] = generate_distribution_band1(n, d, r, sigma, b)
%rng(3) % for reproducibility %i = 5;
X = randn(n,d-1);
X = [ones(n,1) X];
beta = normrnd(0,1,[d,1]);
beta = beta/norm(beta);
beta = beta*b;
%Band Permuted Data
Pi = zeros(n,1);
average = zeros(floor(n/r)+1,d+1);
Y = X*beta + normrnd(0,sigma,[n,1]);
for i = 1:floor(n/r)
   Pi(1 + (i-1)*r: i*r) = randperm(r) + (i-1)*r;
   %blockindex(1+ (i-1)*r: i*r) = i;
   average(i,1:d) = mean(X(1 + (i-1)*r: i*r,:),1);
   average(i,d+1) = mean(Y(1 + (i-1)*r: i*r));
end
average(floor(n/r)+1,1:d) = mean(X(floor(n/r)*r + 1:n,:),1);
average(floor(n/r)+1,d+1) = mean(Y(floor(n/r)*r + 1:n));
average(isnan(average(:,1)),:) = [];
Pi(floor(n/r)*r + 1:n) = randperm(n - floor(n/r)*r) + floor(n/r)*r;
[abc, inv_Pi] = sort(Pi);
Y_B = Y(Pi);
end