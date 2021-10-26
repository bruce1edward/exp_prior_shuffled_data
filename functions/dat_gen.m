% Data generator 
function [X,Y,beta,Pi] = dat_gen(n, d, K, Y_dis, X_dis, b, link)
rng('default') % for reproducibility
%Design for X
if isequal(X_dis,'abs_normal')
  X = abs(randn(n,d));
end
if isequal(X_dis,'normal')
  X = randn(n,d);
end
if isequal(X_dis,'snormal')
  X = pearsrnd(2,1,-1,3,n,d);
end
if isequal(X_dis,'lnormal')
  X = exp(randn(n,d));
end
if isequal(X_dis,'uniform')
  X = rand(n,d);
end
if isequal(X_dis,'chisq')
  X = chi2rnd(1,n,d);
end

%Design for Beta
%b = 1;
beta = abs(normrnd(0,1,[d,1]));
beta = [0;beta]; % intercept
beta = beta/norm(beta)*b;

%Sparse Permuted Data
pi = randperm(K);
Pi = 1:n;
Pi(sort(pi)) = pi;
X_P = [ones(n,1) X(Pi,:)];

%Generate Y
if isequal(Y_dis,'poisson')
    if isequal(link, 'identity')
    mu = X_P*beta;
    Y = poissrnd(mu);
    end
    if isequal(link, 'canonical')
    mu = exp(X_P*beta);
    Y = poissrnd(mu);
    end
    if isequal(link, 'square')
    mu = (X_P*beta).^2;
    Y = poissrnd(mu);
    end
end
if isequal(Y_dis,'binomial')
    p = exp(X_P*beta)./(1 + exp(X_P*beta));
    %Y = binornd(n,p)/n;
    Y = binornd(200000,p);
end
if isequal(Y_dis,'gamma')
    %Assuming the shape parameter is same for yi, to be 1 (exponential case)
    shape = 20;
    scale = exp(X_P*beta);
    Y = gamrnd(shape,scale);
end

if isequal(Y_dis,'normal')
    Y = X_P*beta + randn(n,1)*1;
end

if isequal(Y_dis,'exponential')
    if isequal(link, 'canonical')
     scale = 1./(X_P*beta);
     Y = exprnd(scale);
    else
     scale = exp(X_P*beta);
     Y = exprnd(scale); 
    end
end
end
