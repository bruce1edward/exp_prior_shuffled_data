function [assignments,P] = MAP_est_multi(n, Y, eta, theta, sigma_sq)
Cost = zeros(n);
for i = 1:n
    for j = 1:n
        Cost(i,j) = sum(Y(j,:).*(eta(i,:)))/sigma_sq - (i~=j)*theta;
    end
end
Cost1 = (Cost - min(Cost,[],'all'))*10^6;
[assignments,P] = sparseAssignmentProblemAuctionAlgorithm(Cost1);
end