function [track_para] = EM_DA_mallow_block(Y_permuted, X_permuted, mcmc_steps, m, iter, theta, g_ID)
%X_permuted = X;
uni_g_ID = unique(g_ID);C = cell(numel(uni_g_ID),2);
for g = 1:numel(uni_g_ID)
C{g,1} = find(g_ID == uni_g_ID(g))';C{g,2} = 1:numel(C{g,1});
end
[n,d] = size(X_permuted);
beta_ss = X_permuted\Y_permuted;
sigma_sq_ss = norm(Y_permuted - X_permuted*beta_ss)^2/n;
track_para = zeros(iter, d);
XTX = X_permuted'*X_permuted;
for k = 1:iter
    track = zeros(mcmc_steps/m,n);               % Initialized the hat_Pi\
    track1 = zeros(mcmc_steps/m, d);
    track2 = zeros(mcmc_steps/m, 1);
    Y_hat = X_permuted*beta_ss/sigma_sq_ss; % Computing the estimator of Y_hat
    
    for g = 1:numel(uni_g_ID)
     [track(:,C{g,1}),C{g,2}] = DA(Y_permuted(C{g,1}), Y_hat(C{g,1}), C{g,2}, mcmc_steps, m, theta, C{g,1});
    end
    
    %sample from mixture distribution of beta
    for ii = 1:mcmc_steps/m
        track1(ii,:) = XTX\((X_permuted(track(ii,:),:)')*Y_permuted);
    end
    r = unidrnd(mcmc_steps/m);
    beta_ss = mvnrnd(track1(r,:),inv(XTX)*sigma_sq_ss)';
    
    %sample from mixture distribution of sigma_sq
    r = unidrnd(mcmc_steps/m);
    for ii = 1:mcmc_steps/m
        eta = X_permuted(track(ii,:),:)*beta_ss;
        track2(ii,:) = (Y_permuted - eta)'*(Y_permuted - eta)/(n - d - 1);
    end
    sigma_sq_ss = 1/gamrnd((n - d - 1)/2,2/((n - d - 1)*track2(r)));
    %track_para(k, :) = [beta_ss' sigma_sq_ss];
    track_para(k, :) = beta_ss';
end
end
