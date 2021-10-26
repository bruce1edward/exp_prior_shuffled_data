%MAP
function [assignments] = MAP_est_block_sparse(mu,Y,g_ID,theta)
%mu = Y;theta = 0;
[n,d] = size(Y);
Cost = zeros(n);
blockindex_ = unique(g_ID);
G = length(blockindex_);
index1 = 1:n;
for g = 1:G   
    index = g_ID == blockindex_(g);
    index2 = index1(index);
    n_G = sum(index);
    for i = 1:n_G
        for j = 1:n_G
              Cost(index2(i),index2(j)) = - Y(index2(j))*(mu(index2(i))) - (index2(i)~=index2(j))*theta;
        end
    end
end
min_cost = min(Cost,[],'all');
Cost1 = zeros(n);
for g = 1:G   
    index = g_ID == blockindex_(g);
    index2 = index1(index);
    n_G = sum(index);
    for i = 1:n_G
        for j = 1:n_G
              Cost1(index2(i),index2(j)) = (Cost(index2(i),index2(j)) - min_cost)*10^6;
        end
    end
end

[assignments,P] = sparseAssignmentProblemAuctionAlgorithm(Cost1);
%sum(1:n ~= assignments')
end