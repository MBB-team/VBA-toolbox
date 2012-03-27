function [pe] = ProbError(mus,Qs,ps,in)
% Numerical evaluation of 1D binary classification error
% [pe] = ProbErrorND(mus,Qs,ps,in)
% IN:
%   - mus: cell array of means
%   - qs: cell array of variances
%   - ps: vector of model prior probabilities
% OUT:
%   - pe: the probability of wrongly categorizing a data sample
% SEE ALSO: ProbErrorND, JensenShannon
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

n = length(mus);
try
    in.gri;
catch
    in.gri = -10:1e-3:10;
end
gri = in.gri(:);

p = zeros(length(gri),n);
for i= 1:n
    tmp = mus{i}-gri;
    p(:,i) = exp(-0.5*tmp.^2./Qs{i});
    p(:,i) = ps(i).*p(:,i)./sum(p(:,i));
end
maxp = max(p,[],2);
pe = 1-sum(maxp);
