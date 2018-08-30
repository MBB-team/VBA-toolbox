function [p_BMA] = VBA_BMA (p0, F0)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [p_BMA] = VBA_BMA (p0, F0)
% performs Bayesian Model Averaging (BMA)
%
% IN:
%   - p0: a Kx1 array of VBA posterior structures, which are
%     conditional onto specific generative models
%   - F0: a Kx1 vector of the respective log-model evidences
% OUT:
%   - p_BMA: the resulting posterior structure that describe the marginal 
%     (over models) probability density functions
%
% /////////////////////////////////////////////////////////////////////////

% for retrocompatibility, accept cell array of posterior
if iscell (p0)
    p0 = cell2mat (p0);
end

% shortcuts
% =========================================================================
% number of models
K = length (p0);

% posterior model probabilities
% =========================================================================
ps = VBA_softmax (F0);

% perform averaging
% =========================================================================

% observation parameters
% -------------------------------------------------------------------------
try
    [p_BMA.muPhi, p_BMA.SigmaPhi] = averageMoments ({p0.muPhi}, {p0.SigmaPhi}, ps);
end

% evolution parameters
% -------------------------------------------------------------------------
try
    [p_BMA.muTheta, p_BMA.SigmaTheta] = averageMoments ({p0.muTheta}, {p0.SigmaTheta}, ps);
end

% initial conditions
% -------------------------------------------------------------------------
try
    [p_BMA.muX0, p_BMA.SigmaX0] = averageMoments ({p0.muX0}, {p0.SigmaX0}, ps);
end


% hidden states
% -------------------------------------------------------------------------
try
    % number of timepoints
    T = size (p0(1).muX, 2);
    % initialisation
    mus = cell (K, 1);
    Qs = cell (K, 1);
    % loop over timepoints
    for t = 1 : T
        % collect moments
        for k=1:K
            mus{k} = p0(k).muX(:, t);
            Qs{k} = p0(k).SigmaX.current{t};
        end
        % compute average
        [p_BMA.muX(:, t), p_BMA.SigmaX.current{t}] = averageMoments (mus, Qs, ps);
    end
end

% observation precision
% -------------------------------------------------------------------------
try
    % number of gaussian sources
    nS = numel (p0(1).a_sigma);
    % initialisation
    mus = cell (K, 1);
    Qs = cell (K, 1);
    % loop over sources
    for iS = 1 : nS
        % collect moments
        for k = 1 : K
            mus{k} = p0(k).a_sigma(iS) / p0(k).b_sigma(iS);
            Qs{k}  = p0(k).a_sigma(iS) / p0(k).b_sigma(iS) ^ 2;
        end
        % compute average
        [m, v] = averageMoments (mus, Qs, ps);
        % map to gamma distribution parameters
        p_BMA.b_sigma(iS) = m / v;
        p_BMA.a_sigma(iS) = m * p_BMA.b_sigma(iS);
    end
end

% hidden state precision
% -------------------------------------------------------------------------
% initialisation
mus = cell (K, 1);
Qs = cell (K, 1);
isStochastic = nan (K, 1);
% collect moments, if any 
for k = 1 : K
    try
        mus{k} = p0(k).a_alpha / p0(k).b_alpha;
        Qs{k} = p0(k).a_alpha / p0(k).b_alpha ^ 2;
        isStochastic(k) = ~ isempty(mus{k}) && ~ isinf(mus{k});
    catch
        isStochastic(k) = false;
    end
end
% compute average, if meaninful
% + all deterministic systems
if ~ any (isStochastic) 
    p_BMA.b_alpha = Inf;
    p_BMA.a_alpha = 0;
% + all stochastic systems
elseif all (isStochastic) 
    % average moments
    [m, v] = averageMoments(mus, Qs, ps);
    % map to gamma distribution parameters
    p_BMA.b_alpha = m / v;
    p_BMA.a_alpha = m * p_BMA.b_alpha;
else
% + mixture of stochastic and deterministic systems
    disp('VBA_MBA: Warning! mixture of deterministic and stochastic models!')
    p_BMA.b_alpha = NaN;
    p_BMA.a_alpha = NaN;
end

end



% #########################################################################
% Subfunctions
% #########################################################################

% Compute averages of 1st order (mus) and 2nd order (Qs) moments of
% distributions, weigthed by ps.
function [m, V] = averageMoments (mus, Qs, ps)
% initialisation
V = zeros(size(Qs{1}));
m = zeros(size(mus{1}));
K = length(ps);
% average 1st order moments
for k = 1 : K
    m = m + ps(k) .* mus{k};
end
% average 2nd order moments
for k = 1 : K
    tmp = mus{k} - m;
    V = V + ps(k) .* (tmp * tmp' + Qs{k});
end
end