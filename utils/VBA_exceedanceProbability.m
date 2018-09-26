function ep = VBA_exceedanceProbability (mu, Sigma, options)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% ep = VBA_exceedanceProbability (mu, Sigma, options)
%
% Calculates the exceedance probabilities for mutivariate Gaussian or 
% Dirichlet distributions, i.e. the probability, for each variable, to be
% greater than all the other ones.
%
% IN:
%   - mu/Sigma: sufficient statistics of the pdf
%       -> for a Gaussian distribution, set
%                mu = E[x]
%                Sigma = V[x]
%       -> for a Dirichlet distribution, set 
%               mu = alpha, ie, the Dirichlet pseudo-counts
%               Sigma = [], or NaN, or left undefined
%   - options: structure with the optional fields:
%       + verbose: display textual info {false}
%       + method: 'sampling' or {'analytical'} derivation of the ep
%       + nSamples: number of samples for the sampling method {1e6}
% OUT:
%   - ep: vector of exceedance probabilities
%
% Notes:
% ~~~~~~
%
% The analytical solution for the Dirichlet distribution is based on 
% numerical integration over Gamma distributions [1] derived by 
% Joram Soch [2].
% [1] https://arxiv.org/abs/1611.01439
% [2] mailto:joram.soch@bccn-berlin.de
%
% /////////////////////////////////////////////////////////////////////////

% check parameters
% =========================================================================

% guess distribution
if ~ exist ('Sigma', 'var') || isempty (Sigma) || VBA_isWeird(Sigma)
    form = 'dirichlet';
else
    form = 'gaussian';
    assert (all (size (Sigma) == numel (mu) * [1 1]), ...
            '*** VBA_exceedanceProbability: For Gaussian densities, Sigma must be a n x n matrix if mu has n elements');
end

% fill in defaults
if ~ exist ('options', 'var')
    options = struct;
end
options = VBA_check_struct (options, ...
    'verbose', false, ...
    'method', 'analytical', ...
    'nSamples', 1e6 ...
    );

% catch trivial case to avoid problems with sampling
% =========================================================================
if isempty (mu)
    ep = [];
    return;
end

if isscalar (mu)
    ep = 1;
    return
end

% initialisation
% =========================================================================
K = numel (mu);
ep = ones (K, 1);

% ep computation
% =========================================================================

switch form
    % ---------------------------------------------------------------------
    case 'gaussian'
        
        switch options.method
            case 'sampling'
                r_samp = VBA_random ('Gaussian', mu, Sigma, options.nSamples);
                [~, j] = max (r_samp);
                tmp = histcounts (j, 0.5 : K + 0.5);
                ep = tmp' / options.nSamples;
        
            case 'analytical'
                c = [1; -1];
                for k = 1 : K
                    for l = setdiff (1 : K, k)
                        ind = [k, l];
                        m = VBA_vec(mu(ind));
                        V = Sigma(ind, ind);
                        ep(k) = ep(k) * VBA_PPM (c' * m, c' * V *c, 0);
                    end
                end
                ep = ep ./ sum (ep);
        end
        
    case 'dirichlet'
    % ---------------------------------------------------------------------
        switch options.method
            case 'sampling'
                r_samp = VBA_random ('Dirichlet', mu, options.nSamples);
                [y, j]  = max(r_samp);
                % remove failed samples in limit cases
                if any (isnan (VBA_vec (y))) 
                    j(isnan (y)) = []; 
                    options.nSamples = numel (j);
                    warning ('VBA_exceedanceProbability: unstable parametrization, only %d%% of samples were correctly generated.', round(100*Nsamp/numel(y)));
                end
                tmp = histc (j, 1 : length (mu));
                ep = tmp / options.nSamples;
                
            case 'analytical'
                for k = 1 : K
                    f = @(x) integrand (x, mu(k), mu(1 : K ~= k));
                    ep(k) = integral (f, eps, Inf);
                end
                ep = ep' ./ sum (ep);
        end
     
end

end

% #########################################################################
% Integrand function for numerical integration of the Dirichlet ep
function p = integrand (x, aj, ak)
    p = ones(size (x));
    % Gamma CDF
    for k = 1 : numel (ak)
        p = p .* gammainc (x, ak(k));                         
    end
    % Gamma PDF
    p = p .* exp ((aj - 1) .* log (x) - x - gammaln (aj));         
end               

