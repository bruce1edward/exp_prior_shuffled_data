function [Betahat, Xi, conv] = blockcoord(X, Y, lambda, tol, maxit)

if  nargin < 5
    maxit = 10000;
end
if nargin < 4
   tol = 1E-4; 
end

conv = zeros(1, maxit);

[n, m] = size(Y);
[Q, T] = qr(X, 0);
%[Betacur, Residual] = R \ Y;
Yhat = Q * (Q' * Y);
R = Y - Yhat;
R_rownorms = sqrt(sum(R.^2, 2));
Xi = zeros(n, m);
it = 1;
conv(it) = sum(R_rownorms.^2)/(n * m);

while true
   
   R =  Y - Yhat;
   
   R_rownorms = sqrt(sum(R.^2, 2));
   shrinkage = 1 - (n * m * lambda ./ (2 * sqrt(n) * R_rownorms));
   shrinkage = shrinkage .* (shrinkage > 0);
    
   %for i=1:n
   Xi = R .* repmat(shrinkage, [1 m]) / sqrt(n); 
   %end
   
   
   Z = Y - Xi * sqrt(n);

   Yhat = Q * (Q' *  Z);
   
   it = it + 1;
   
   conv(it) = sum(sum((Z - Yhat).^2))/(n * m) + lambda * sum(sqrt(sum(Xi.^2, 2)));
   
   if conv(it-1) - conv(it) < tol
       break;
   end
   
   if it >= maxit
      break; 
   end
      
end

conv = conv(1:it);


Betahat = T \ Q'*(Y - sqrt(n) * Xi);

end



