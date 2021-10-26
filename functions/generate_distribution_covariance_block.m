% Data generator covariance
function [x,y,Pi,sigma,blockindex,Q] = generate_distribution_covariance_block(n, p, q, sigmax, sigmay, rho, k, r)
%x = randn(n,p)*sqrt(sigma2x);
%S = cov(x);
%y = randn(n,q)*sqrt(sigma2y);
mu = zeros(p+q,1);
%sigma = [eye(p)*sigmax ones(p,q)*rho; ones(q,p)*rho eye(q)*sigmay];
sigma = ones(p+q,p+q)*rho;
for i = 1:p+q
    if i<=p
       sigma(i,i) = sigmax;
    else
       sigma(i,i) = sigmay; 
    end
end
s = mvnrnd(mu,sigma,n);
x = s(:,1:p);
y = s(:,p+1:p+q);
blockindex = zeros(n,1);
Q = zeros(n);
for i = 1:n/r
   p_i = randperm(r, k);
   P_i = 1:r;
   P_i(sort(p_i)) = p_i;
   Pi(1 + (i-1)*r: i*r) = P_i + (i-1)*r;
   blockindex(1+ (i-1)*r: i*r) = i;
   Q(1 + (i-1)*r: i*r,1 + (i-1)*r: i*r) = (1 - k/r - (k/r)/(r-1)) * eye(r) + ((k/r)/(r-1)) * ones(r);
end
end

