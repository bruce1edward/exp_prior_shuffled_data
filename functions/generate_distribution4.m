function [X_c,Y_s,Pi,inv_Pi,blockindex] = generate_distribution4(n, d, r, noise, b)
%rng(3) % for reproducibility
X = randn(n,d);
%beta = normrnd(0,1,[d,1]);
%beta = beta/norm(beta);
%beta = beta*b;
Y = sin(X - pi); %+ normrnd(0,noise,[n,1]);
[Y_s, inv_s] = sort(Y);
X_c = X(inv_s);
%block Permuted Data
Pi = zeros(n,1);
blockindex = zeros(n,1);
for i = 1:n/r
   Pi(1 + (i-1)*r: i*r) = randperm(r) + (i-1)*r;
   blockindex(1+ (i-1)*r: i*r) = i;
end
[abc, inv_Pi] = sort(Pi);
end