function X1 = sampleDist(f,M,N,b)
% SAMPLEDIST  Sample from an arbitrary distribution
%     sampleDist(f,M,N,b) retruns an array of size X of random values 
%     sampled from the distribution defined by the probability density 
%     function refered to by handle f, over the range b = [min, max].  
%     M is the threshold value for the proposal distribution, such that 
%     f(x) < M for all x in b.
%  
%     sampleDist(...,true) also generates a histogram of the results
%     with an overlay of the true pdf.
%  
%     Examples: 
%     %Sample from a step function over [0,1]:
%     X = sampleDist(@(x)1.3*(x>=0&x<0.7)+0.3*(x>=0.7&x<=1),...
%                    1.3,1e6,[0,1],true);
%     %Sample from a normal distribution over [-5,5]:
%     X = sampleDist(@(x) 1/sqrt(2*pi) *exp(-x.^2/2),...
%                    1/sqrt(2*pi),1e6,[-5,5],true);
%
n1 = 0;
%N = 1000;b = [0,5*log(n)];M = 5; 
X1 = zeros(N,1);
counter = 0;

while n1 < N && counter < 1000
    x = b(1) + rand(5*N,1)*diff(b);
    uM = M*rand(5*N,1);
    %x = b(1) + rand(1000,1)*diff(b);
    %uM = M*rand(1000,1);
    x = x(uM < f(x));
    %length(x)
    if isempty(x)
        error('No Points Sampled')
    end
    X1(n1+1:min([n1+length(x),N])) = x(1:min([length(x),N - n1]));
    n1 = n1 + length(x);
    counter = counter+1;
end
end