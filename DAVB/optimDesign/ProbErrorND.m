function [pe] = ProbErrorND(mus,Qs,ps,in)
% Numerical evaluation of ND binary classification error
% [pe] = ProbErrorND(mus,Qs,ps,in)
% IN:
%   - mus: cell array of ND-means
%   - qs: cell array of NDxND covariance matrices
%   - ps: vector of model prior probabilities
% OUT:
%   - pe: the probability of wrongly categorizing a data sample
% SEE ALSO: ProbError, JensenShannon
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

nm = length(mus);
n = length(mus{1});
try
    in.gri;
catch
    in.gri = -10:2e-1:10;
end
gri = in.gri(:);

% get vectors of all ND-grid points
[C] = get_nkdraws(length(gri),n,1);

% Evaluate Gausian densities on the ND-grid
np = size(C,2);
p = zeros(np,nm);
for i= 1:nm
    iV = pinv(Qs{i});
    mu =  mus{i};
    for j=1:np
        tmp = mu-gri(C(:,j));
        p(j,i) = exp(-0.5*tmp'*iV*tmp);
    end
    p(:,i) = ps(i).*p(:,i)./sum(p(:,i));
end

% Derive selection error probability
maxp = max(p,[],2);
pe = 1-sum(maxp);
