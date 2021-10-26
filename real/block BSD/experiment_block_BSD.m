function [] = experiment_block_BSD(which)
rng(which)
addpath('../../functions') 
k__n = 0.05:0.05:0.4;
k_n = k__n(which);
data1 = readtable('../../data/BSD/data.csv');
X = data1(:,2:12);
X = X{:,:};
X(:,8) = round(X(:,8)*47 - 8);
data = readtable('../../data/BSD/data1.csv');
Y = data(:,2);
X1 = data(:,3:14);
Y = Y{:,:};
X1 = X1{:,:};
X1 = [X1 X1(:,3).*X1(:,6:8)];
[n,d] = size(X1);
iter = 100;
iters = 100;
mcmc_steps = 2000;
burn_steps = 500;

g_ID = findgroups(X(:,5)); %0.7514/167
error_est = zeros(iters,5);
dev_est = zeros(iters,5);

for i = 1 : iters
% Data generator block
P = data_generator_block1(n, g_ID, k_n);
K_n = sum(1:n~=P)/n;
Y_P = Y(P);
dev = @(beta) 2*nansum(Y.*log(Y./exp([ones(n,1) X1]*beta)) - (Y - exp([ones(n,1) X1]*beta)));
%Oracle
glm1 = fitglm(X1,Y,'linear','Distribution','poisson');
beta_oracle = glm1.Coefficients.Estimate;
%Naive
glm2 = fitglm(X1,Y_P,'linear','Distribution','poisson');
beta_naive = glm2.Coefficients.Estimate;
%-------------------------------
%EM
beta_EM = EM_mal_Block_GLM(Y_P, X1, iter, mcmc_steps, burn_steps, beta_naive, g_ID);
%EMM
beta_EMM = EM_mal_Block_GLM1(Y_P, X1, iter, mcmc_steps, burn_steps, beta_naive, g_ID, Y);
%------------------------------
uni_g_ID = unique(g_ID); Q_g = zeros(n,n);
for g = 1:numel(uni_g_ID)
index_g = find(g_ID == uni_g_ID(g));
Q_g(index_g,index_g) = ones(numel(index_g),numel(index_g))*(k_n/(n-1));
for gg = 1:numel(index_g)
    Q_g(index_g(gg),index_g(gg)) = 1 - k_n;
end
end
%LL
[beta_LL, norm_score_LL, converged_LL] = LL_Armijo(n,d,X1,Q_g,Y_P,beta_naive);
%Chamber 
[beta_C, norm_score_C, converged_C] = Chamber_Armijo(n,d,X1,Q_g,Y_P,beta_naive);

error_est(i,:) = [norm(beta_naive - beta_oracle)/norm(beta_oracle), norm(beta_EM  - beta_oracle)/norm(beta_oracle), norm(beta_EMM  - beta_oracle)/norm(beta_oracle), norm(beta_LL  - beta_oracle)/norm(beta_oracle), norm(beta_C  - beta_oracle)/norm(beta_oracle)];
dev_est(i,:) = [dev(beta_naive) dev(beta_EM) dev(beta_EMM) dev(beta_LL) dev(beta_C)];
end
save(['results_experiment_bsd_' num2str(which) '.mat'],'error_est','dev_est','K_n');
exit
end