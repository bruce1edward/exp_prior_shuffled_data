function [hat_Pi,order,count] = mcmc_covariance1(x, y, order, mcmc_steps, burn_steps, omega_xy, theta)
%{
n = 1000; p = 3; q = 3; sigmax = 1; sigmay = 1; rho = 0.5;k = floor(0.2*n); 
[x,y,Pi] = generate_distribution_covariance(n, p, q, sigmax, sigmay, rho, k);
z_naive = [x(Pi,:), y];S_naive = cov(z_naive)*(n-1)/n;
order = 1:n; mcmc_steps = 8000; burn_steps = 4000;iters = 200;omega = S_naive\eye(p+q);omega_xy = omega(1:p,p+1:p+q);
theta = log(n-k);
%}
    [n,p] = size(x);
    [n,q] = size(y);
    %Pi_x = zeros(n,p);
    hat_Pi = zeros(n,n);
    %Pi_y = zeros(n,q);
    %track = zeros(mcmc_steps,2);
    count = 0;
    for m_step = 1: mcmc_steps
        %sum(order ~=1:n)
                order_ = order;        % replace the permutation with the old one
                i = randi(n);          % generate random position 1
                j = randi(n);          % generate random position 2
                order_(i) = order(j);  % swap the position 1 with 2
                order_(j) = order(i);  % swap the position 2 with 1
                x_diff = repmat((x(order(i),:) - x(order(j),:))', [1,q]);
                y1 = repmat(y(i,:),[p,1]);
                y2 = repmat(y(j,:),[p,1]);
                %A1 = x_diff.*y1 -  x_diff.*y2;
                A = x_diff.*y1 -  x_diff.*y2;
                %check
                %y'*(x(order,:) - x(order_,:))/n
                %(x(order,:)' - x(order_,:)')*y/n
                %S_pi = [zeros(p,p), A1; A1' zeros(q,q)]; % compute the acceptance
                %p22 = trace(omega*S_pi)/2 + theta*(sum(order ~= 1:n) - sum(order_ ~= 1:n)); % compute the acceptance raito
                %p21 = trace(omega_xy*A') + theta*(sum(order ~= 1:n) - sum(order_ ~= 1:n)); % compute the acceptance raito
                p21 = sum(omega_xy.*A,'all') + theta*(sum(order ~= 1:n) - sum(order_ ~= 1:n));
                %track(m_step,:) = [p21 p22];
                 if 0 < p21              % if acceptance raito is > 0, then accpet
                       order = order_;
                       count = count + 1;
                 else
                    if rand < exp(p21)   % if not then don't 
                       order = order_; 
                       count = count + 1;
                    end
                 end
                if m_step > burn_steps   % the burning steps
                    %order_t = zeros(n,1);
                    %order_t(order) = 1:n;
                    %Pi_x = Pi_x + x(order,:)/(mcmc_steps - burn_steps);
                    %Pi_y = Pi_y + y(order_t,:)/(mcmc_steps - burn_steps);
                    lin_index = (1:n) + n*(order-1);  % linear indexing
                    hat_Pi(lin_index) = hat_Pi(lin_index) + 1;  %update the estimating Pi
                end
    end
end