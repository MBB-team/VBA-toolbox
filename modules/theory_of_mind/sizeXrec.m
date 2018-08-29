function nx= sizeXrec(i,P)
% returns the number of k-ToM's hidden sates
% function nx = sizeXrec(i,P)
% Marie Devaine wrote this in November 2015 (comments: JD).
% The formula follows from an analysis of the recursive sequence.
% IN:
%   - i: k-ToM's sophistication level
%   - P: total nb of evol/obs params
% OUT:
%   - nx: number of k-ToM's hidden sates

if i==0
    nx = 2;
else
    i_1 = i-1;
    if i_1>0
        power2 = exp(((i_1-1):-1:0).*log(2));
        S = 2^i + power2*((1:i_1)*(3*P+2)-1)';
    else
        S = 2;
    end
    nx = S+ i*(3*P+2)-1;
end
end