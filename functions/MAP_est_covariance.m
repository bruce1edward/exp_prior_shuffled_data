function [assignments] = MAP_est_covariance(n, omega_xy ,x, y, theta)
%x = x_observed;y = y_center;theta = 0;
xx = -x*omega_xy*y';Cost = zeros(n);
for i = 1:n
    for j = 1:n
        Cost(i,j) = xx(i,j) - (i~=j)*theta;
    end
end
Cost1 = (Cost - min(Cost,[],'all'))*10^6;
[assignments,P] = sparseAssignmentProblemAuctionAlgorithm(Cost1);
end