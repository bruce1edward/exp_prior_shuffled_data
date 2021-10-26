%Chamber_Armijo
function [betacur, norm_score, converged] = Chamber_Armijo(n,d,X,Q,ystar,beta_start)
%beta_start = beta_naive;ystar = Y_permuted;
X = [ones(n,1) X];
maxiter = 10000;
norm_res = zeros(maxiter,1);
iter = 0;
tol = 10^-6;
tau = 0.05;
gamma = 0.8;
flag = false;
betacur = beta_start;

score = @(beta) X'*(ystar - Q*exp(X*beta));
hess = @(beta) - X'*(Q*repmat(exp(X*beta), [1 (d+1)]).* X);
foo = @(s) sqrt(sum(s.^2));

rhs = score(betacur);
iter = iter + 1;
norm_res(iter) = foo(rhs);
converged = 1*(norm_res(iter) < tol);
track = zeros(maxiter,1);

while(converged~= 1 && iter < maxiter)
    score_xcur = score(betacur);
    hess_xcur = hess(betacur);
    normsq_score_xcur = foo(score_xcur);
    track(iter) = cond(hess(betacur));
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
    if  flag
        break; 
    end
    betacur = betacur - gamma^m * (hess_xcur\score_xcur);
    iter = iter + 1;
    rhs = score(betacur);
    norm_res(iter) = foo(rhs);
    converged = 1*(norm_res(iter) < tol);
end
norm_score = norm_res(iter);
%plot(norm_res,'*')

end