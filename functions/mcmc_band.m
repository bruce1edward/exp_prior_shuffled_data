function [hat_Pi,order] = mcmc_band(Y, Y_hat, order, n, mcmc_steps, burn_steps, theta, C)
    hat_Pi = zeros(n,n);               % Initialized the hat_Pi
    %Pi_Y = zeros(n,1);
    for m_step = 1: mcmc_steps
                order_ = order;        % replace the permutation with the old one
                i = randi(n);          % generate random position 1
                j = randi(n);          % generate random position 2
                order_(i) = order(j);  % swap the position 1 with 2
                order_(j) = order(i);  % swap the position 2 with 1
                A = Y_hat(order_,:) - Y_hat(order,:); % compute the acceptance 
                p21 = sum(Y(:) .* A(:)) + theta*(trace(C(order,:)) - trace(C(order_,:))); % compute the acceptance raito 
                 if 0 < p21              % if acceptance raito is > 0, then accpet
                       order = order_;
                 else
                    if rand < exp(p21)   % if not then don't 
                       order = order_; 
                    end
                 end
                if m_step > burn_steps   % the burning steps
                    %order_t = zeros(n,1);
                    %order_t(order) = 1:n;
                    %Pi_Y = Pi_Y + Y(order_t)/(mcmc_steps - burn_steps);
                    lin_index = (1:n) + n*(order-1);  % linear indexing
                    hat_Pi(lin_index) = hat_Pi(lin_index) + 1;  %update the estimating Pi
                end
    end
end

