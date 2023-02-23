n = 10; d = 1; 
X =  (1:1:10)'; b = 1; sigma = 0.1;
rng(1)
Y = X*b + sigma*randn(n,1);
%plot(X,Y,'*');
Y_P = Y;
Y_P(3) = Y(8);
Y_P(8) = Y(3);
plot(X,Y_P,'*')
P_true = [1 2 8 4 5 6 7 3 9 10];
%----------------
Pi = 1:10;
norm_density = @(Pi) exp(-sum((Y - X(Pi)*b).^2)/(2*sigma^2));
norm_density_shuffled = @(Pi) exp(-sum((Y_P - X(Pi)*b).^2)/(2*sigma^2));

%check
%psi__(0.1)
%poisscdf(n,exp(0.1) - 1)*exp(exp(0.1) - 1)

gamma_density = @(gamma, Pi) exp(- gamma*sum(Pi ~= 1:n))/(factorial(10)*exp(-gamma*10)*fun_psi__(gamma));
gamma_density(0.1, 1:n)

%Exact Expectation
P = perms(1:n);
np = factorial(n);
p_density = zeros(np,1);
for j = 1 : np
p_density(j) = norm_density_shuffled(P(j,:)) * gamma_density(5, P(j,:));
end
%p_density =  p_density - max(p_density);
%p_density = softmax(p_density);
P((p_density > 0.1),:)
[A,B] = sort(p_density);

top11_P = P(B((np - 10):np),:);
for i = 1:11
   hamming_dis(i) = sum(top11_P(i,:)~=P_true);
end

plot(hamming_dis, "*")

%----------fixing Pi maximizing--------------
Pi = 1 : n;
f_p_y = @(gamma) norm_density_shuffled(Pi) * gamma_density(gamma, Pi);
gamma = fminbnd(f_p_y,0,1000);
%------------------------
Pi1 = P_true;
f_p_y1 = @(gamma) norm_density_shuffled(Pi1) * gamma_density(gamma, Pi1);
n_f_p_y1 = @(gamma) (norm_density_shuffled(Pi1) * gamma_density(gamma, Pi1))*(-1);
gamma1 = fminbnd(n_f_p_y1,0,10);

grid = 0:0.1:10;
for i = 1:101
    fun_val(i) = f_p_y1(grid(i));
end
plot(grid, fun_val,'*')
%--------------------------

%Exact Expectation
sum_psi1 = @(gamma) fun_sum_psi(gamma, norm_density_shuffled, gamma_density);
gamma1 = fminbnd(sum_psi1,0,10);

grid = 0:0.1:10;
for i = 1:101
    fun_val(i) = sum_psi1(grid(i));
end
sum_psi1(3)

plot(grid(1:(i-1)), fun_val(1:(i-1)),'*')

[A,B] = min(fun_val(1:(i-1)));
abc = grid(1:(i - 1));
abc(B)
