%Sorting 
function [pi] = E_Pi2(eta,Y)
[n,d] = size(Y);
[A,index] = sort(Y,'descend');
[B,index1] = sort(eta);
index_inv(index) = 1:n;
pi = index1(index_inv);
end
