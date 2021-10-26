%C1
function [betacur] = Chamber(n,d,X,Q,ystar,beta_start)
%Q = Q2;beta_start = beta_naive;ystar = Y;X = X1;
maxiter = 10000;
norm_res = zeros(maxiter,1);
iter = 0;
tol = 10^-6;
X = [ones(n,1) X];
X_ = Q*X;
betacur = beta_start;
mucur = exp(X*betacur);  % start from the naive solution
rescur = ystar - Q*mucur;     % current residual 
iter = iter + 1;
rhs = X' * rescur;
norm_res(iter) = sqrt(sum(rhs.^2));
converged = 1*(norm_res(iter) < tol);

while(converged~= 1 && iter < maxiter)
    % compute quantities for the linear system
    w =  mucur; % weights (variances)
    %plot(w,'*')
    WX = repmat(w, [1 d+1]).* X;
    betacur = X_'*WX\(X'*rescur + X_'*WX * betacur);
    %quadratic problem (reverse normalize equation)
    mucur = exp(X * betacur);  % start from the naive solution
    rescur = ystar - Q*mucur;     % current residual 
    iter = iter + 1;
    rhs = X'* rescur;
    norm_res(iter) = sqrt(sum(rhs.^2));
    converged = 1*(norm_res(iter) < tol);
end

%plot(norm_res,'*')

end