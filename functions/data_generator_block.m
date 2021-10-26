function [P] = data_generator_block(n,g_ID,k_n)
blockindex_ = unique(g_ID);
G = length(blockindex_);
P = 1:n;
%index2 = 1:n;
%index1 = [];
%track1 = zeros(G,1);%max(track1)
for g = 1:G   
    %g = 11;
    index = g_ID == blockindex_(g);
    n_G = sum(index);
    %track1(g) = n_G;
    if n_G < 10
        P_ = randperm(n_G);
        P__ = P(index);
        P(index) = P__(P_);
        %sum(P ~= 1:n)
    end
    P_ = randperm(n_G, floor(k_n*n_G));
    P__ = P(index);
    P__(sort(P_)) = P__(P_);
    P(index) = P__;
end
%sum(track1(track1 > 10))/n
end