function [psi_gamma] = psi__(gamma)
psi_gamma = 0;
for i = 0:10
    psi_gamma = psi_gamma + (exp(gamma) - 1)^i/factorial(i);
end
end