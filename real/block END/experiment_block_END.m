function [] = experiment_block_END(seed)
%seed = 7;
rng(seed)
addpath('../../functions') 
k__n = 0.1:0.05:0.4;
iteration = 1:100;
[gr1, gr2] =  ndgrid(1:numel(k__n), 1:numel(iteration));
nposs = numel(k__n) * numel(iteration);
Gr = [reshape(gr1, [nposs 1]), reshape(gr2, [nposs 1])];
Gr_row = Gr(seed,:);
k_n = k__n(Gr_row(1));
%Data(Blocking)
data = readtable('../../data/END/elnino_l.csv','ReadVariableNames',false);
Y = data(2:end,11);X = data(2:end,[8:10 12]);lat_lon = data(2:end,[6 7]);
Y = Y{:,:};X = X{:,:};lat_lon = lat_lon{:,:};
Y = cellfun(@str2num,Y);X = cellfun(@str2num,X);lat_lon = cellfun(@str2num,lat_lon);
[n,d] = size(X);[n,m] = size(Y);X = [ones(n,1) X];d = d + 1;
X_wo = X(:,2:end);X_centered = X_wo - repmat(mean(X_wo,1), [n 1]);
g_ID = findgroups(lat_lon(:,1),lat_lon(:,2));%numel(unique(g_ID))/numel(Y)
% Data generator block
P = data_generator_block(numel(Y), g_ID, k_n);%sum(P~=1:n)/n
K_n = sum(P ~= 1:numel(Y))/numel(Y);
Y_permuted = Y(P);
%Naive
beta_naive = X\Y_permuted;
%Oracle
beta_oracle = X\Y;
%centered resposne
Y_permuted1 = Y_permuted - mean(Y_permuted);
%EM Block
iter = 200; mcmc_steps = 2000; burn_steps = 500;
[beta_EM_block, sigma_sq_EM_block] = EM_mal_Block(Y_permuted, X, iter, mcmc_steps, burn_steps, 0, beta_naive, g_ID);
%EM 
iter = 200; mcmc_steps = 2000; burn_steps = 500; 
[beta_EM_block1, sigma_sq_EM_block1] = EM_mal_Block1(Y_permuted, X, iter, mcmc_steps, burn_steps, beta_naive, g_ID, Y);
%LL
uni_g_ID = unique(g_ID);X_centered_Q = zeros(n,d-1);
for g = 1:numel(uni_g_ID)
index_g = find(g_ID == uni_g_ID(g));
len = numel(index_g);
if len == 1
    Q_g = 1;
else
    Q_g = ones(len,len)*(K_n/(n-1));
    for gg = 1:len
       Q_g(gg,gg) = 1 - K_n;
    end
end
X_centered_Q(index_g,:) = Q_g*X_centered(index_g,:);
end
beta_LL = X_centered_Q\Y_permuted1;
beta_LL = [mean(Y_permuted - X_wo*beta_LL) ; beta_LL];
%Chamber
beta_C = (X_centered'*X_centered_Q)\(X_centered'*Y_permuted1);
beta_C = [mean(Y_permuted - X_wo*beta_C) ; beta_C];

error_est = [norm(beta_naive - beta_oracle)/norm(beta_oracle), norm(beta_EM_block  - beta_oracle)/norm(beta_oracle), norm(beta_EM_block1  - beta_oracle)/norm(beta_oracle), norm(beta_LL  - beta_oracle)/norm(beta_oracle), norm(beta_C  - beta_oracle)/norm(beta_oracle)];
r_square = 1 - [norm(Y - X*beta_naive)^2/norm(Y- mean(Y))^2, norm(Y - X*beta_EM_block)^2/norm(Y- mean(Y))^2 ,norm(Y - X*beta_EM_block1)^2/norm(Y- mean(Y))^2, norm(Y - X*beta_LL)^2/norm(Y- mean(Y))^2, norm(Y - X*beta_C)^2/norm(Y- mean(Y))^2];
save(['results_experiment_END_' num2str(Gr_row(1) + (Gr_row(2)-1)*7) '.mat'],'error_est','r_square');
exit
end
