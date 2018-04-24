function [DJS,b,muy,Vy] = VBA_JensenShannon(mus,Qs,sources,w,base)
% evaluates the Jensen-Shannon divergence (DJS)
% function [DJS,b] = VBA_JensenShannon(mus,Qs,sources,w,base)
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
%   - sources: structure indicating distribution type
%   - w: nx1 vector of weights
%   - base: base for the log mapping ({'2'}, '10', or 'e')
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

if ~exist('sources','var')
    sources = struct('type', 0); % default to gaussian
end
sources = VBA_check_struct(sources, ...
            'out', 1:numel(mus{1}) ...
);

if ~exist('w','var')
    w = ones(numel(mus),1);
    w = w/sum(w);
end

if ~exist('base','var')
    base = '2';
end

% _________________________________________________________________________
% catch unimplemented cases 

if numel(sources) > 1 % TODO: implement approximation for multisource models
    error('*** Jensen-Shannon divergence cannot be computed (yet) for multisource models');
end

if sources.type == 2 % TODO: compute JS divergence for multinomial models
    error('*** Jensen-Shannon divergence cannot be computed (yet)for multinomial models');
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
    
% _________________________________________________________________________
% compute the divergence


% get mixture of entropy 
%
% $$ \sum_{i=1}^n p_i H(P_i) $$
%

sH = 0;
for i=1:n
  % get weighted sum of entropies 
  switch sources.type
      
    % gaussian
    case 0 
      [e] = eig(full(Qs{i}));
      logDet = sum(log_b(e));
      sH = sH + 0.5*w(i).*logDet;
      
    %binomial
    case 1 
      sH = sH -sum(mus{i}.*log_b(mus{i})) -sum((1-mus{i}).*log_b(1-mus{i}));
      
  end
end

% get entropy of mixture:
%
% $$ H(\sum_{i=1}^n p_i P_i) $$
%

Hy = 0;
Vy = zeros(p);

% Compute moments of the mixture distribution
% - 1st order moment
muy = horzcat(mus{:}) * w ;

switch sources.type
    
 % gaussian
 case 0
    % get second order moment of sum of densities
    for i=1:n
        tmp = mus{i} - muy;
        tmp = tmp*tmp' + Qs{i};
        Vy = Vy + w(i).*tmp;
    end
    % get Gaussian approx entropy
    [e] = eig(full(Vy));
    Hy = 0.5*sum(log_b(e));
    
  % binomial 
  case 1
    Hy = -sum(muy.*log_b(muy)) -sum((1-muy).*log_b(1-muy));
    
end

% get Jensen-Shannon approximation
DJS = Hy - sH;

% _________________________________________________________________________
% get error probability upper bound
Hp = -sum(w.*log_b(w));
b = max([-Inf,Hp - DJS]);

end


    
    
    