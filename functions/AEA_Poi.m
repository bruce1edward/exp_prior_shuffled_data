%Alternative Estimating Algorithm (AEA)
function [beta_hat, xi_hat] = AEA_Poi(X,Y,lambda,beta_start)
%Debug
%profile on 
%Y = Y_P;
%lambda = lambda1;beta_start = beta_naive;
[n,d] = size(X);
xi_hat = zeros(n,1);
epsilon = 10^-5;
%maxiter = 1000;
fval = zeros(10000,1);
X1 = [ones(n,1) X];
%Loss = @(coef, f) -dot(X*coef + f,Y) + sum(exp(X*coef + f)) + lambda*sum(abs(f));
%coef = beta_start;xi = 0;
Loss = @(coef, xi) 1/n*(-dot(X1*coef + sqrt(n)*xi,Y) + sum(exp(X1*coef + sqrt(n)*xi))) + lambda*sum(abs(xi));
iter = 2;
fval(iter) = Loss(beta_start,0);
fval(iter-1) = 0;
beta_hat = beta_start;
while norm(fval(iter) - fval(iter-1)) > epsilon %&& iter < maxiter
    beta_cur = beta_hat;
    eta = X1*beta_cur;
  %STH update
  STH = @(eta) 1.*(Y < exp(eta) - sqrt(n)*lambda).*(log(Y + sqrt(n)*lambda) - eta) + 1.*(Y > exp(eta) + sqrt(n)*lambda).*(log(Y - sqrt(n)*lambda) - eta);
  xi_hat = STH(eta)/sqrt(n);
  %norm(xi_hat1 - xi_hat);
  %{
   for i = 1 : n  
        if Y(i) < exp(eta(i)) - lambda
         xi_hat(i) = log(Y(i)+lambda) - eta(i);
        elseif Y(i) <=  exp(eta(i)) + lambda
         xi_hat(i) = 0;
        else 
         xi_hat(i) = log(Y(i) - lambda) - eta(i);
        end  
    end
  %}
  
  %One step IRLS update
    XWX = X1' * (repmat(exp(eta + sqrt(n)*xi_hat), [1 d + 1]) .* X1);
    beta_hat = beta_cur + XWX\(X1'*(Y - exp(eta + sqrt(n)*xi_hat)));
      %line search
      %Warm start
      %multiple steps IRLS (1 > 3)
    iter = iter + 1;
    fval(iter) = Loss(beta_hat,xi_hat);
end
%norm(beta_hat - beta_true)
%norm(beta_hat - beta_naive)
%iter
%plot(xi_hat,'*')
%plot(fval(1:10))
%plot(diff(fval(1:6)),'*')
%profile off
%profile viewer
end