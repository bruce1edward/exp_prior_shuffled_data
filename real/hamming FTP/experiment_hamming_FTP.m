function [] = experiment_hamming_FTP(which)
%which = 1;
rng(which)
addpath('..\..\functions')
load('..\..\data\FTP\ftp_processed.mat')
K_n = 0.1 : 0.05 : 0.4;
iteration = 1:100;
[gr1, gr2] =  ndgrid(1:numel(K_n), 1:numel(iteration));
nposs = numel(K_n) * numel(iteration);
Gr = [reshape(gr1, [nposs 1]), reshape(gr2, [nposs 1])];
Gr_row = Gr(which,:);k_n = K_n(Gr_row(1));
iter = 200; % number of iteration for EM alogrithm
mcmc_steps = 4000; burn_steps = 2000; 
error_est = zeros(1,4);
%------------------------------------
%Sparse Permtuation 
K = floor(n*k_n); 
pi = randperm(K);
Pi = 1:n;
Pi(sort(pi)) = pi;
%Naive 
z_naive = [X(Pi,:) , Y];S_naive = cov(z_naive)*(n-1)/n;R_naive = corrcoef(z_naive);
%Oracle
z_oracle = [X , Y];S_oracle = cov(z_oracle)*(n-1)/n;R_oracle = corrcoef(z_oracle);
%Centered x and y
X_center = X - mean(X,1);
Y_center = Y - mean(Y,1);
%Variance Matrix
var_vector = sqrt(diag(S_oracle)); sigmaxy = zeros(d+m,d+m);
for ii = 1:d+m
    for jj = 1:d+m
        sigmaxy(ii,jj) = var_vector(ii)*var_vector(jj);
    end
end
%EM
order = 1:n;
omega_start = S_naive\eye(d+m); 
S_pi = EM_covariance(X_center(Pi,:), Y_center, order, mcmc_steps, burn_steps, iter, omega_start, d, m, 0);
%EMM
order = 1:n;
omega_start = S_naive\eye(d+m); 
theta = Choose_theta(n,K);
S_pi_mallows = EM_covariance(X_center(Pi,:), Y_center, order, mcmc_steps, burn_steps, iter, omega_start, d, m, theta);
%Robust
sig_robust = robustcov(z_naive);

error_est(1) = norm(R_naive - R_oracle, 'fro');
error_est(2) = norm(S_pi./sigmaxy - R_oracle, 'fro');
error_est(3) = norm(S_pi_mallows./sigmaxy - R_oracle, 'fro');
error_est(4) = norm(sig_robust./sigmaxy - R_oracle, 'fro');
save(['results_experiment_FTP_' num2str(Gr_row(1) + (Gr_row(2)-1)*7) '.mat'],'error_est');
exit
end

