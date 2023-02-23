rng('default')
addpath('functions') 
addpath('data')
load('Italian_survey_data_linkage.mat');X_link = X;X = X(:,1);Y_link = Y;Y = Y(:,1);
n = numel(X);X = [ones(n,1) X];K_n = 0.4;
%sum(Pi~=1:n)
pi_ = randperm(floor(K_n*numel(Y)));Pi1 = 1:numel(Y);Pi1(sort(pi_)) = pi_;
Y_permuted = Y(Pi1);
%Naive
beta_naive = X\Y_permuted;
sigma_sq_naive = norm(Y_permuted - X*beta_naive)^2/n;
%Oracle
beta_oracle = X\Y;
sigma_sq_oracle = norm(Y - X*beta_oracle)^2/n;
%DA
iter = 1000;mcmc_steps = 4000; % mcmc steps
m = 40;order = 1:n;
track_para = DA_mallow(Y_permuted, X, order, mcmc_steps, mcmc_steps/m, iter, 0);
%DA regulaized
iter = 1000;mcmc_steps = 4000; % mcmc steps
m = 40;order = 1:n;theta = Choose_theta(n,floor(K_n*n));
track_para1 = DA_mallow(Y_permuted, X, order, mcmc_steps, mcmc_steps/m, iter, theta);
%DA HB
iter = 1000;mcmc_steps = 4000; % mcmc steps
m = 40;order = 1:n;
track_para2 = DA_HB(Y_permuted, X, order, mcmc_steps, mcmc_steps/m, iter, log(n), 1);

iter = 200; % number of iteration for EM alogrithm
mcmc_steps = 8000; burn_steps = 4000; 
%EM
order = 1:n;
[beta_EM, sigma_EM] = EM_mal_tricks(Y_permuted, X, iter, mcmc_steps, burn_steps, 0, beta_naive, order);
%EMM
order = 1:n;
theta = Choose_theta(n,floor(K_n*n));
[beta_EMM, sigma_EMM] =  EM_mal_tricks(Y_permuted, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);
%EM Empirical Bayes
order = 1:n; theta = log(n);
[beta_EMMB, sigma_EMMB] =  EM_mal_EB(Y_permuted, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);


[f1,xi1] = ksdensity(track_para(:,1)); 
[f11,xi11] = ksdensity(track_para1(:,1)); 
[f111,xi111] = ksdensity(track_para2(:,1)); 

figure
hold on
p1 = plot(xi1,f1,'r','LineWidth', 1.5);
p2 = plot(xi11,f11,'k','LineWidth', 1.5);
p21 = plot(xi111,f111,'g','LineWidth', 1.5);
p3 = xline(beta_oracle(1),'b','LineWidth', 1.5);
p4 = xline(beta_naive(1),'y','LineWidth', 1.5);
p5 = xline(beta_EMMB(1),'m','LineWidth', 1.5);
p6 = xline(beta_EMM(1),'c','LineWidth', 1.5);

xlabel('\beta_0', 'FontSize', 14)
ylabel('Posterior distribution', 'FontSize', 14)
lh = legend([p1 p2 p21 p3 p4 p5 p6], {'Without', 'With', 'HB','Oracle','Naive', 'EB', 'EMH'}, 'location', 'best');
lh.NumColumns = 1;
lh.FontSize = 10;
lh.ItemTokenSize = [30 70];  
legend('boxoff')
grid on
set(gcf, 'Color', 'w');
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 7;
hold off

export_fig('DA_example1.pdf')


[f2,xi2] = ksdensity(track_para(:,2),'kernel','box'); 
[f22,xi22] = ksdensity(track_para1(:,2),'kernel','box'); 
[f222,xi222] = ksdensity(track_para2(:,2),'kernel','box'); 
figure
hold on
p1 = plot(xi2,f2,'r','LineWidth', 1.5);
p2 = plot(xi22,f22,'k','LineWidth', 1.5);
p21 = plot(xi222,f222,'g','LineWidth', 1.5);
p3 = xline(beta_oracle(2),'b','LineWidth', 1.5);
p4 = xline(beta_naive(2),'y','LineWidth', 1.5);
p5 = xline(beta_EMMB(2),'m','LineWidth', 1.5);
p6 = xline(beta_EMM(2),'c','LineWidth', 1.5);
xlabel('\beta_1', 'FontSize', 14)
ylabel('Posterior distribution', 'FontSize', 14)
lh = legend([p1 p2 p21 p3 p4 p5 p6], {'Without', 'With', 'HB','Oracle','Naive', 'EB', 'EMH'}, 'location', 'best');
lh.NumColumns = 1;
lh.FontSize = 10;
lh.ItemTokenSize = [30 70];  
legend('boxoff')
grid on
set(gcf, 'Color', 'w');
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 7;
hold off

export_fig('DA_example2.pdf')

[f3,xi3] = ksdensity(track_para(:,3)); 
[f33,xi33] = ksdensity(track_para1(:,3)); 
[f333,xi333] = ksdensity(track_para2(:,3)); 
figure
hold on
p1 = plot(xi3,f3,'r','LineWidth', 1.5);
p2 = plot(xi33,f33,'k','LineWidth', 1.5);
p21 = plot(xi333,f333,'g','LineWidth', 1.5);
p3 = xline(sigma_sq_oracle,'b','LineWidth', 1.5);
p4 = xline(sigma_sq_naive,'y','LineWidth', 1.5);
p5 = xline(sigma_EMMB,'m','LineWidth', 1.5);
p6 = xline(sigma_EMM,'c','LineWidth', 1.5);
xlabel('\sigma^2', 'FontSize', 14)
ylabel('Posterior distribution', 'FontSize', 14)
lh = legend([p1 p2 p21 p3 p4 p5 p6], {'Without', 'With', 'HB','Oracle','Naive', 'EB', 'EMH'}, 'location', 'best');
lh.NumColumns = 1;
lh.FontSize = 10;
lh.ItemTokenSize = [30 70];  
legend('boxoff')
grid on
set(gcf, 'Color', 'w');
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 7;
hold off

export_fig('DA_example3.pdf')