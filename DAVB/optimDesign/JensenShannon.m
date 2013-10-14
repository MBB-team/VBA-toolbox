function [DJS,b,muy,Vy] = JensenShannon(mus,Qs,ps,binomial,base)
% evaluates the Jensen-Shannon divergence (DJS)
% function [DJS,b] = JensenShannon(mus,Qs,ps,binomial)
% This function evaluates the DJS:
% - either from a set of N-D Gaussian densities, in which case those are
% defined through their first- and second-order moments,
% - or binomial densities, in which case the first-order moments are
% sufficient.
% In addition, a set of individual weights for each of the n component
% densities should be provided.
% IN:
%   - mus: nx1 cell array of 1st-order moments
%   - Qs: nx1 cell array of 2nd-order moments
%   - ps: nx1 vector of weights
%   - binomial: flag for binomial densities {0}
%   - base: base for the log mapping (can be '2' or '10') {'10'}
% OUT:
%   - DJS: the Jensen-Shannon divergence
%   - b: the associated lower-bound on the ensuing probability of
%   classicafication error 
%   - muy/Vy: the 1st- and 2d- order moments of the Laplace approx to the
%   mixture
% SEE ALSO: ProbError, ProbErrorND

try,binomial;catch,binomial=0;end
try,base;catch,base='2';end

n = length(mus); % number of models to be compared
if ~binomial
    Vy = zeros(size(Qs{1}));
end
muy = zeros(size(mus{1}));
sH = 0;

for i=1:n
    % get first order moment of mixture
    muy = muy + ps(i).*mus{i};
    % get weighted sum of entropies 
    if binomial
        sH = sH -sum(mus{i}.*LOG(mus{i},base)) -sum((1-mus{i}).*LOG((1-mus{i}),base));
    else
        [e] = eig(full(Qs{i}));
        logDet = sum(LOG(e,base));
        sH = sH + 0.5*ps(i).*logDet;
    end
end

% get mixture entropy
if binomial
    Hy = -sum(muy.*LOG(muy,base)) -sum((1-muy).*LOG((1-muy),base));
else
    % get second order moment of sum of densities
    for i=1:n
        tmp = mus{i} - muy;
        tmp = tmp*tmp' + Qs{i};
        Vy = Vy + ps(i).*tmp;
    end
    % get Gaussian approx entropy
    [e] = eig(full(Vy));
    Hy = 0.5*sum(LOG(e,base));
end

% get Jensen-Shannon approximation
DJS = Hy - sH;

% get error probability upper bound
Hp = -sum(ps.*LOG(ps,base));
b = max([-Inf,Hp - DJS]);


function lo = LOG(x,base)
switch base
    case '2'
        lo = log2(x);
    case '10'
        lo = log(x);
end

    
    
    
    