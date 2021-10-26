function [Q] = generate_Q1(n, blockindex, k_n, r)
% blockindex = g_ID(idx(1:round(0.8*n)));
blockindex_ = unique(blockindex);
G = length(unique(blockindex));
B = zeros(G,n);
Q = zeros(n);
index_ = 1:n;
for g = 1:n/r       % iterate over group 
   Q = (1 - k_n - k_n/(r-1)) * eye(r) + (k_n/(r-1)) * ones(r);
end
end