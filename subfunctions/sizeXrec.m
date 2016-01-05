function fx= sizeXrec(i,P)
% sizeRec returns the numbers of hidden sates corresponding to a particular
% ToM-level (i) associated wit the number (P) of parameters of the generative
% model (volatility, temperature, bias...)
% The formula follows from an analysis of the recursive sequence.
% Marie Devaine wrote this in November 2015
%
if i==0
    fx=2;
else
    i_1=i-1;
    if i_1>0
        power2= exp(((i_1-1):-1:0).*log(2));
        S=2^i + power2*((1:i_1)*(3*P+2)-1)';
    else
        S=2;
    end
    fx=S+ i*(3*P+2)-1;
end
end