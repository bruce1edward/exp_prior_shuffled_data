function [track_para] = DA_mallow(Y_permuted, X_permuted, order, mcmc_steps, m, iter, theta)
[n,d] = size(X_permuted);
beta_ss = X_permuted\Y_permuted;
sigma_sq_ss = norm(Y_permuted - X_permuted*beta_ss)^2/n;
track_para = zeros(iter, d+1);
XTX = X_permuted'*X_permuted;
for k = 1:iter
    track = zeros(mcmc_steps/m,n);               % Initialized the hat_Pi\
    track1 = zeros(mcmc_steps/m, d);
    track2 = zeros(mcmc_steps/m, 1);
    count = 0;
    Y_hat = X_permuted*beta_ss/sigma_sq_ss; % Computing the estimator of Y_hat
    for m_step = 1: mcmc_steps
                order_ = order;        % replace the permutation with the old one
                i = randi(n);          % generate random position 1
                j = randi(n);          % generate random position 2
                order_(i) = order(j);  % swap the position 1 with 2
                order_(j) = order(i);  % swap the position 2 with 1
                A = Y_hat(order_,:) - Y_hat(order,:); % compute the acceptance 
                p21 = sum(Y_permuted(:) .* A(:)) + theta*(sum(order ~= 1:n) - sum(order_ ~= 1:n)); % compute the acceptance raito 
                 if 0 < p21              % if acceptance raito is > 0, then accpet
                       order = order_;
                 else
                    if rand < exp(p21)   % if not then don't 
                       order = order_; 
                    end
                 end
                if mod(m_step,m) == 0    % the burning steps
                    count = count + 1;
                    track(count,:) = order;
                end
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
    track_para(k, :) = [beta_ss' sigma_sq_ss];
    %track_para(k, :) = beta_ss';
end
end
