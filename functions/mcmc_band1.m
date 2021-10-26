function [hat_Pi,order,effective_mcmc_step,effective_count] = mcmc_band1(Y, Y_hat, order, n, mcmc_steps, burn_steps, r)
    hat_Pi = zeros(n,n);               % Initialized the hat_Pi
    effective_mcmc_step = 0;
    effective_count = 0;
    for m_step = 1: mcmc_steps
                order_ = order;        % replace the permutation with the old one
                i = randi(n);          % generate random position 1
                j = randi([max(i-r,1), min(i+r,n)],1,1); % generate random position 2
                if abs(order(i) - order(j)) > r % Check the band constraint
                    effective_mcmc_step = effective_mcmc_step + 1;
                    continue
                end
                order_(i) = order(j);  % swap the position 1 with 2
                order_(j) = order(i);  % swap the position 2 with 1
                A = Y_hat(order_,:) - Y_hat(order,:); % compute the acceptance 
                p21 = sum(Y(:) .* A(:)); % compute the acceptance raito 
                 if 0 < p21              % if acceptance raito is > 0, then accpet
                       order = order_;
                 else
                    if rand < exp(p21)   % if not then don't 
                       order = order_; 
                    end
                 end
                if m_step > burn_steps - effective_mcmc_step   % the burning steps
                    lin_index = (1:n) + n*(order-1);  % linear indexing
                    hat_Pi(lin_index) = hat_Pi(lin_index) + 1;  %update the estimating Pi
                    effective_count = effective_count + 1;
                end
    end
end

