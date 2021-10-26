function [X,Y_B,beta,Pi,average] = generate_distribution_band_GLM1(n, d, r, b)
%rng(3) % for reproducibility
X = randn(n,d);
X_w_intercept = [ones(n,1) X];
%X = rand(n,d-1)*2*sqrt(3) - sqrt(3);
beta = normrnd(0,1,[d + 1,1]);
beta = beta/norm(beta);
beta = beta*b;
%Band Permuted Data
Pi = zeros(n,1);
average = zeros(floor(n/r)+1,d+1);
for i = 1:floor(n/r)
   Pi(1 + (i-1)*r: i*r) = randperm(r) + (i-1)*r;
   average(i,1:d) = mean(X(1 + (i-1)*r: i*r,:),1);
end
Pi(floor(n/r)*r + 1:n) = randperm(n - floor(n/r)*r) + floor(n/r)*r;
average(floor(n/r)+1,1:d) = mean(X(floor(n/r)*r + 1:n,:),1);
mu = exp(X_w_intercept(Pi,:)*beta);
Y_B  = poissrnd(mu);
for i = 1:floor(n/r)
   average(i,d+1) = mean(Y_B(1 + (i-1)*r: i*r));
end
average(floor(n/r)+1,d+1) = mean(Y_B(floor(n/r)*r + 1:n));
end