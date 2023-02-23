function [p_density] = sum_psi(gamma, norm_density_shuffled, gamma_density)
n = 10;
P = perms(1:n);
np = factorial(n);
p_density = 0;
for j = 1 : np
p_density = p_density + norm_density_shuffled(P(j,:)) * gamma_density(gamma, P(j,:));
end
p_density = -log(p_density);
end