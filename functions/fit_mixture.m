function [betahat, alphahat, sigmahat, conv] = fit_mixture(X, y, control)
%X = X_centered;y = Y_permuted1; control = "LS";
%b = betahat; alph = alphahat; sigma = sigmahat;
    tol = 1E-4;
    maxit = 100;
    [n,p] = size(X);
    tausq = mean(y.^2);
    if control == "LS"
        betahat = X\y;
        sigmahat = sqrt(sum((y - X*betahat).^2)/(n-p));
        alphahat = 0.5; 
    end
    
    if control == "robust"
        [betahat,stat] = robustfit(X(:,2:end),y,'huber');
        sigmahat = stat.mad_s;
        %betahat = beta_robust(2);
        %sigmahat = stat.mad_s;
        betahat_LS = X\y;
        alphahat = max(min(.99, sqrt(sum(betahat_LS.^2))/sqrt(sum(betahat.^2))),0.01);
        %alphahat = 0.5;
    end

    %%% computation of pseudo-EM estimator
    % negative log-likelihood
    nloglik = @(b, alph, sigma) mean(-log(alph * (1/sqrt(2 * pi * sigma^2)) * exp(-(y - X*b).^2/(2 * sigma^2))...
        + (1-alph) * (1/sqrt(2 * pi * tausq)) * exp(-y.^2 / (2 * tausq))));

    it = 1;
    conv = zeros(maxit+1, 1);        
    nloglik_cur = nloglik(betahat, alphahat ,sigmahat);
    conv(it) = nloglik_cur;

    while it < maxit
        r = X * betahat - y;
        
        pnum =  alphahat * (1/sqrt(2 * pi * sigmahat^2)) * exp(-r.^2/(2*sigmahat^2));
        pdenom = pnum + (1 - alphahat) * (1/sqrt(2 * pi * tausq)) * exp(-y.^2 / (2 * tausq));
        pcur = pnum ./ pdenom;
        w = pcur;
       
        Xw = X.*repmat(sqrt(w),1,p); 
        yw = sqrt(w).* y;

        betahatnew = Xw\yw;        
        sigmahatnew = sqrt(sum((yw - Xw * betahatnew).^2) / sum(w));
        alphahatnew = mean(pcur);
        nloglik_new = nloglik(betahatnew, alphahatnew, sigmahatnew);

        if nloglik_cur - nloglik_new < -tol^2
            break
        else
            betahat = betahatnew;
            alphahat = alphahatnew;
            sigmahat = sigmahatnew;
            it = it+1;
            conv(it) = nloglik_new;
            if nloglik_cur - nloglik_new < tol
                break
            else
                nloglik_cur = nloglik_new;
            end
        end
    end
  end