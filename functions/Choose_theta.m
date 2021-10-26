function [theta] = Choose_theta(n,K)
%n = numel(C{g,1}); K = K_N; K = floor(k_n*n)
%n = r ; K = k;
upbound1 = @(theta) poisscdf(n-K,exp(theta))/poisscdf(n,exp(theta)-1) - 0.001;
if n == K
    theta = 0;
else
    theta = fzero(upbound1,round(log(n - K),2));
    if isnan(theta)
        theta = log(n-K);
    end
end
%{
theta = 0 : 0.01 : 10;
%theta = 6;
num = zeros(numel(theta),1);
for i = 1 : numel(theta)
    num(i) = upbound1(theta(i));
end
figure
hold on
p1 = plot(theta, num,'k*');
xlabel('\theta')
ylabel('Prob')
title([' n = ' num2str(n) ' k = ' num2str(floor(K))] , 'Interpreter','Latex')
%legend([p1 p2], {'Approximate','Exact'}, 'location', 'best');
hold off
%}
end