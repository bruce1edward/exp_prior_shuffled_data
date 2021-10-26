rng('default')
addpath('..\..\functions')
r_n = 3:1:10;nmethod = 5;
iters = 100;iter = 400; % number of iteration for EM alogrithm
mcmc_steps = 10000; burn_steps = 2000; 
error_est = zeros(numel(r_n),iters,nmethod);
r_square = zeros(numel(r_n),iters,nmethod);

for i = 1:numel(r_n)
for j = 1:iters
n = 1000; d= 20; r = floor(r_n(i)); sigma = 1; b = 3;
[X,Y_B,beta,Pi,average] = generate_distribution_band_GLM1(n, d, r, b);
%Naive 
glm = fitglm(X,Y_B,'linear','Distribution','poisson'); 
beta_naive = glm.Coefficients.Estimate;
%Oracle
glm1 = fitglm(X(Pi,:),Y_B,'linear','Distribution','poisson'); 
beta_oracle = glm1.Coefficients.Estimate;
%EM
order = 1:n;
beta_EM = EM_mal_GLM(n, Y_B, X, iter, mcmc_steps, burn_steps, 0, beta_naive, order);
%EMM
order = 1:n;
beta_EMM = EM_mal_GLM_band(n, Y_B, X, iter, mcmc_steps, burn_steps, r, beta_naive, order);
%average
glm2 = fitglm(average(:,1:end-1),average(:,end),'linear','Distribution','poisson'); 
beta_average = glm2.Coefficients.Estimate;

error_est(i,j,1) = norm(beta_naive - beta)/b;
error_est(i,j,2) = norm(beta_oracle  - beta)/b;
error_est(i,j,3) = norm(beta_EM  - beta)/b;
error_est(i,j,4) = norm(beta_EMM  - beta)/b;
error_est(i,j,5) = norm(beta_average  - beta)/b;
end
end

save('Local_GLM.mat');
