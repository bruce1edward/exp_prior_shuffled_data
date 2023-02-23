function [hat_Pi] = exact_EM(Y_P, X, b, sigma, n, gamma)
%gamma = log(n);
%----------------
norm_density_shuffled = @(Pi) (1/sqrt(2*sigma^2))*exp(-sum((Y_P - X(Pi)*b).^2)/(2*sigma^2));
gamma_density = @(gamma, Pi) exp(- gamma*sum(Pi ~= 1:n))/(factorial(n)*exp(-gamma*n)*poisscdf(n,exp(gamma) - 1)*exp(exp(gamma) - 1));
%check
hat_Pi = zeros(n,n);               % Initialized the hat_Pi
P = perms(1:n);
np = factorial(n);
p_density = zeros(np,1);
for j = 1 : np
p_density(j) = norm_density_shuffled(P(j,:)) * gamma_density(gamma, P(j,:));
lin_index = (1:n) + n*(P(j,:)-1);  % linear indexing
hat_Pi(lin_index) = hat_Pi(lin_index) + 1*p_density(j);  %update the estimating Pi
end
hat_Pi = hat_Pi/sum(p_density);
end

