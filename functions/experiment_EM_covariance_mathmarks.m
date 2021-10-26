function [] = experiment_EM_covariance_mathmarks(which)
rng(which)%which = 1;
k_n = [0.2, 0.4, 0.6];
mcmc_steps = 500:500:8000; 
iter = 20;
[gr1, gr2] =  ndgrid(1:numel(k_n), 1:numel(mcmc_steps));
nposs = numel(mcmc_steps)*numel(k_n);
Gr = [reshape(gr1, [nposs 1]), reshape(gr2, [nposs 1])];
Gr_row = Gr(which,:);

data = importdata('..\data\mathmarks\mathmarks.dat');
[n,d] = size(data.data);
p = 2;x = data.data(:,1:p);q = d - p;y = data.data(:,p+1:d);

mc_steps = mcmc_steps(Gr_row(2));
bu_steps = mc_steps*0.25;

K = floor(k_n(Gr_row(1))*n);
theta = [0 1./(20:-1:1)*Choose_theta(n,K) Choose_theta(n,K):0.4:Choose_theta(n,K)+2];
gridtheta = numel(theta);
error_est = zeros(iter,gridtheta);
error_est_naive = zeros(iter,1);
for i = 1:iter
%Sparse Permtuation 
pi = randperm(K);Pi = 1:n;Pi(sort(pi)) = pi;
%Naive 
z_naive = [x(Pi,:) , y];S_naive = cov(z_naive)*(n-1)/n; R_naive = corrcoef(z_naive);
%Oracle
z_oracle = [x , y];S_oracle = cov(z_oracle)*(n-1)/n; R_oracle = corrcoef(z_oracle);
error_est_naive(i) = norm(R_naive - R_oracle, 'fro');
%Variance Matrix
var_vector = sqrt(diag(S_oracle)); sigmaxy = zeros(d,d);
for ii = 1:d
    for jj = 1:d
        sigmaxy(ii,jj) = var_vector(ii)*var_vector(jj);
    end
end
%Centered x and y
x_center = x - mean(x,1);
y_center = y - mean(y,1);
for j = 1:gridtheta
order = 1:n;
omega_start = S_naive\eye(d); 
S_pi_mallows = EM_covariance(x_center(Pi,:), y_center, order, mc_steps, bu_steps, 200, omega_start, p, q, theta(j));
error_est(i,j) = norm(S_pi_mallows./sigmaxy - R_oracle, 'fro');
end
end

save(['C:\Users\zwang39\OneDrive - George Mason University\Paper(EM)\code\nips2021\functions\diag\result\results_experiment_' num2str(which) '.mat']);

MSE_beta = mean(error_est,1);
std_beta = std(error_est,1);

figure('visible', 'off');
hold on 
p1 = errorbar(theta, MSE_beta, 3*std_beta/sqrt(iter), 'k-*', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'k');
xlabel('theta', 'FontSize', 14)
xlim([min(theta) max(theta)])
ax=gca;
ax.FontSize = 14;
xticks(0:1:Choose_theta(n,K) + 2)
%set(gca,'xticklabel',[])
p2 = yline(mean(error_est_naive), 'LineWidth', 1.5);
title([' k/n = ' num2str(k_n(Gr_row(1))) ' , ' ' mcmc steps = ' num2str(mc_steps)])
lh = legend([p1 p2], {'EMM(\theta)','Naive'}, 'location', 'best');
grid on
set(gcf, 'Color', 'w');
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 7;
hold off

export_fig(['diag/plot/results_experiment_' num2str(which) '.pdf'])

end