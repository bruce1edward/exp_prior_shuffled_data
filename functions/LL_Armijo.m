%LL_Armijo
function [betacur, norm_score, converged] = LL_Armijo(n,d,X,Q,ystar,beta_start)
%beta_start = beta_naive;ystar = Y;
X = [ones(n,1) X];
maxiter = 10000; 
norm_res = zeros(maxiter,1);
iter = 0;
tol = 10^-6;
X_ = Q*X;
betacur = beta_start;
rescur = ystar - Q*exp(X*betacur);     % current residual 
rhs = X_' * rescur;
iter = iter + 1;
norm_res(iter) = sqrt(sum(rhs.^2));
converged = 1*(norm_res(iter) < tol);
tau = 0.05;
gamma = 0.3;
flag = false;
score = @(beta) X_'*(ystar - Q*exp(X*beta));
hess = @(beta) - X_'*(Q*repmat(exp(X*beta), [1 (d+1)]).* X);
foo = @(s) sqrt(sum(s.^2));
%track = zeros(maxiter,1);
%track1 = zeros(maxiter,1);
while(converged~= 1 && iter < maxiter)
    score_xcur = score(betacur);
    hess_xcur = hess(betacur);
    normsq_score_xcur = foo(score_xcur);
    %track(iter) = cond(hess(betacur));
    % check Armijo rule
    m = 0;
    while foo(score_xcur) - foo(score(betacur - gamma^m * (hess_xcur\score_xcur))) < tau * gamma^m * normsq_score_xcur 
    %while foo(score_xcur) - foo(score(betacur - gamma^m * (hess_xcur\score_xcur))) < tau * normsq_score_xcur 
        m = m+1;
        if gamma^m < tol^2
            flag = true;
           break;
        end
    end
    if flag 
        break; 
    end
    betacur = betacur - gamma^m * (hess_xcur\score_xcur);
    iter = iter + 1;
    rhs = score(betacur);
    norm_res(iter) = foo(rhs);
    converged = 1*(norm_res(iter) < tol);
end
norm_score = norm_res(iter);
%plot(norm_res, '*')
end