function [y,labels] = generateBMM(n,alpha,lambda)
%
% function y = generateGMM(n,lambda,mu)
% This function generates n observations following a MoB distribution
% given by the class frequencies alpha, and their average pattern lambda


[D,K] = size(lambda);
labels = zeros(n,K);

k = 1;
sumAlphak = alpha(1);
% seed = 1e4;
for i = 1:n
    ratio = i./n;
    if ratio>=sumAlphak && i ~= n
        k = k+1;
        sumAlphak = sumAlphak+alpha(k);
    end
    for d=1:D
        [y(d,i)] = VBA_sampleFromArbitraryP([lambda(d,k),1-lambda(d,k)]',[1,0]',1);
    end
    labels(i,k) = 1;
end
