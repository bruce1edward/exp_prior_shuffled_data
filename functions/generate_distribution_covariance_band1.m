% Data generator covariance
function [x,y,Pi,sigma,average] = generate_distribution_covariance_band1(n, p, q, sigmax, sigmay, rho, r)
%x = randn(n,p)*sqrt(sigma2x);
%S = cov(x);rho = 0.8;
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
Pi = zeros(n,1);
average = zeros(floor(n/r)+1,p+q);
for i = 1:floor(n/r)
   Pi(1 + (i-1)*r: i*r) = randperm(r) + (i-1)*r;
   average(i,1:p) = mean(x(1 + (i-1)*r: i*r,:),1);
   average(i,p+1:p+q) = mean(y(1 + (i-1)*r: i*r,:),1);
end
Pi(floor(n/r)*r + 1:n) = randperm(n - floor(n/r)*r) + floor(n/r)*r;
average(floor(n/r)+1,1:p) = mean(x(floor(n/r)*r + 1:n,:),1);
average(floor(n/r)+1,p+1:p+q) = mean(y(floor(n/r)*r + 1:n,:),1);
average(isnan(average(:,1)),:) = [];
end

