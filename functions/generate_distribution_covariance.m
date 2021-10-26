% Data generator covariance
function [x,y,Pi,sigma] = generate_distribution_covariance(n, p, q, sigmax, sigmay, rho, k)
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
pi = randperm(k);
Pi = 1:n;
Pi(sort(pi)) = pi;
end

