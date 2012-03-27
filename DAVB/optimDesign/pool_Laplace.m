function [mu,S,pf] = pool_Laplace(m,Q,families,pm)
% derive prior predictive densities over families of models
% function [mu,S] = pool_Laplace(m,Q,families)
% IN:
%   - m: nmX1 cell array of prior means
%   - Q: nmX1 cell array of prior covariance matrices
%   - families: nfX1 cell aray of indices of models belonging to each
%   family or partition of model space
%   - pm: nmX1 vector of models prior probabilities
% OUT:
%   - mu: nfX1 cell array of prior means
%   - S: nfX1 cell array of prior covariance matrices
%   - pf: nfX1 vector of families prior probabilities
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

nf = length(families);
nm = length(m);
IND = zeros(nm,1);
for i=1:nm
    for j=1:nf
        if ismember(i,families{j})
            IND(i) = j;
        end
    end
end
if ~exist('pm','var') || isempty(pm)
    pm = zeros(nm,1);
    for i=1:nm
        n = length(find(IND==IND(i)));
        pm(i) = 1./n;
    end
    pm = pm./sum(pm);
end

% 1- derive means
mu = cell(nf,1);
pf = zeros(nf,1);
for i=1:nf
    in = families{i};
    pi = pm(in);
    pf(i) = sum(pi);
    pi = pi./sum(pi);
    mu{i} = 0;
    for j=1:length(in)
        mu{i} = mu{i} + pi(j).*m{in(j)};
    end
end

% 2- derive covariance matrices
S = cell(nf,1);
for i=1:nf
    in = families{i};
    pi = pm(in);
    pi = pi./sum(pi);
    S{i} = 0;
    for j=1:length(in)
        tmp = m{in(j)} - mu{i};
        S{i} = S{i} + pi(j).*( tmp*tmp' + Q{in(j)} );
    end
end


