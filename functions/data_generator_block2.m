function [P] = data_generator_block2(n, Y,g_ID,k_n)
gridtheta = 500;
gamma = linspace(0,1000,gridtheta);
K_N = zeros(gridtheta,1);
for i = 1:gridtheta
P_map = MAP_est_block_sparse(Y,Y,g_ID,gamma(i));
K_N(i) = sum(P_map' ~= 1:n)/n;
end
[abc, pos] = min(abs(K_N - k_n));
P = MAP_est_block_sparse(Y,Y,g_ID,gamma(pos));
end