% Data generator block
function [X,Y,Y_P,beta,Pi,blockindex,Q] = generate_distribution_block(n, d, r, k, sigma, b)
X = randn(n,d-1);
X = [ones(n,1) X];
beta = normrnd(0,1,[d,1]);
beta = beta/norm(beta);
beta = beta*b;
Y = X*beta + normrnd(0,sigma,[n,1]);
Pi = zeros(n,1);
blockindex = zeros(n,1);
Q = zeros(n);
%block Permuted Data
for i = 1:n/r
   p_i = randperm(r,k);
   P_i = 1:r;
   P_i(sort(p_i)) = p_i;
   Pi(1 + (i-1)*r: i*r) = P_i + (i-1)*r;
   blockindex(1+ (i-1)*r: i*r) = i;
   Q(1 + (i-1)*r: i*r,1 + (i-1)*r: i*r) = (1 - k/r - (k/r)/(r-1)) * eye(r) + ((k/r)/(r-1)) * ones(r);
end
%{
p_i = randperm(n,floor((k/r)*n));
P_i = 1:n;
P_i(sort(p_i)) = p_i;
Y_P = Y(P_i);
Pi = P_i;
Q = (1 - (k/r) - (k/r)/(n-1)) * eye(n) + ((k/r)/(n-1)) * ones(n);
%[abc, inv_Pi] = sort(Pi);
%}
Y_P = Y(Pi);
end

