function [DJS,b,muy,Vy] = VBA_JensenShannon(mus,Qs,w,binomial,base)
% evaluates the Jensen-Shannon divergence (DJS)
% function [DJS,b] = VBA_JensenShannon(mus,Qs,ps,binomial)
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
%   - w: nx1 vector of weights
%   - binomial: flag for binomial densities {0}
%   - base: base for the log mapping (can be '2' or '10') default = 2
% OUT:
%   - DJS: the Jensen-Shannon divergence
%   - b: the associated lower-bound on the ensuing probability of
%   classicafication error 
%   - muy/Vy: the 1st- and 2d- order moments of the Laplace approx to the
%   mixture
%
% /////////////////////////////////////////////////////////////////////////
%
% Given a set of probability distributions P_1, ..., P_n and weights $p$,
% the Jensen-Shannon Divergence JSD is:
%
% $$ 
%   JSD_p(P_1, ..., P_n) 
%     = H(\sum_{i=1}^n w_i P_i) - \sum_{i=1}^n w_i H(P_i) 
%     = Hy - sH
% $$
%
% where $H(P)=E[-log(P)]$ is the Shannon entropy for the distribution P 
%
% /////////////////////////////////////////////////////////////////////////

% _________________________________________________________________________
% check inputs
try,binomial;catch,binomial=0;end

if ~exist('base','var')
    base = '2';
end

% _________________________________________________________________________
% initialization

% short hand for log computation in the desired base
switch base
    case '2'
        log_b = @(x) log2(x);
    case '10'
        log_b = @(x) log10(x);
    case 'e'
        log_b = @(x) log(x);
end

% number of models to be compared
n = length(mus); 
% data size
p = numel(mus{1});

% Compute moments of the mixture distribution
% - 1st order moment
muy = horzcat(mus{:}) * w ;
% - 2nd order moment
    
% _________________________________________________________________________
% compute the divergence


% get mixture of entropy 
%
% $$ \sum_{i=1}^n p_i H(P_i) $$
%
% For each density, as the sources are independents, the entropy of the 
% joint distribution over sources is equal to the sum of source-wise 
% entropies: 
% $$ H[P(y)] = H[P(y_s1)] + H[P(y_s2)] $$
%
sH = 0;
for i=1:n
    % get weighted sum of entropies 
    if binomial
        sH = sH -sum(mus{i}.*log_b(mus{i})) -sum((1-mus{i}).*log_b(1-mus{i}));
    else
        [e] = eig(full(Qs{i}));
        logDet = sum(log_b(e));
        sH = sH + 0.5*w(i).*logDet;
    end
end

% get entropy of mixture
%
% $$ H(\sum_{i=1}^n p_i P_i) $$
%
% as sources are independent,
Hy = 0;
if binomial
    Hy = -sum(muy.*log_b(muy)) -sum((1-muy).*log_b(1-muy));
else
    % get second order moment of sum of densities
    for i=1:n
        tmp = mus{i} - muy;
        tmp = tmp*tmp' + Qs{i};
        Vy = Vy + w(i).*tmp;
    end
    % get Gaussian approx entropy
    [e] = eig(full(Vy));
    Hy = 0.5*sum(log_b(e));
end




% get Jensen-Shannon approximation
DJS = Hy - sH;

% _________________________________________________________________________
% get error probability upper bound
Hp = -sum(w.*log_b(w));
b = max([-Inf,Hp - DJS]);


end


    
    
    