function ep = VBA_exceedanceProbability (distribName, varargin)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% ep = VBA_exceedanceProbability ('Gaussian', mu, Sigma, [options])
% ep = VBA_exceedanceProbability ('Dirichlet', alpha, [options])
%
% Calculates the exceedance probabilities for mutivariate Gaussian or 
% Dirichlet distributions, i.e. the probability, for each variable, to be
% greater than all the other ones.
%
% IN:
%   - distribName: type of distribution. Could be 'Gaussian' or 'Dirichlet'
%
%   - mu, Sigma, or alpha: sufficient statistics of the pdf
%       -> for a Gaussian distribution
%                mu = E[x]
%                Sigma = V[x]
%       -> for a Dirichlet distribution 
%               alpha = the Dirichlet pseudo-counts
%
%   - options: optional structure with the optional fields:
%       + verbose: display textual info {false}
%       + method: 'sampling' or {'analytical'} derivation of the ep
%       + nSamples: number of samples for the sampling method {5e6}
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

% check distribution and parse moments
switch distribName
    case 'Gaussian'
        mu = varargin{1};
        Sigma = varargin{2};
        K = numel (mu);
        assert (isvector (mu) && ismatrix (Sigma) && all (size (Sigma) == numel (mu) * [1 1]), ...
            '*** VBA_exceedanceProbability: For Gaussian densities, mu should be a n vector and Sigma a n x n matrix.');
    case 'Dirichlet'
        alpha = varargin{1};
        K = numel (alpha);
        assert (isvector(alpha) && all (alpha > 0), ...
            '*** VBA_exceedanceProbability: For Dirichlet densities, alpha should be a positive vector.');
    otherwise
        error('*** VBA_exceedanceProbability: unknown probability density type.');
end

% check options
parser = inputParser;
parser.PartialMatching = false;

parser.addParameter ('verbose', false, @islogical);
parser.addParameter ('method', 'analytical', @(z) ismember(z, {'analytical','sampling'}));
parser.addParameter ('nSamples', 1e7, @(z) z > 0 && mod(z,1) == 0);

if isstruct(varargin{end})
    parser.parse (varargin{end});
else
    parser.parse ();
end
options = parser.Results ;

% catch trivial cases
% =========================================================================
if K < 2
    ep = NaN;
    return;
end

% ep computation
% =========================================================================

switch distribName
    % ---------------------------------------------------------------------
    case 'Gaussian'
        
        switch options.method
            
            case 'sampling'
                r_samp = VBA_random ('Gaussian', mu, Sigma, options.nSamples);
                ep = sum(bsxfun (@eq, r_samp, max (r_samp)), 2);
                ep = ep / sum (ep);
                
            case 'analytical'
                if ~ VBA_isdiag (Sigma)
                    error('*** VBA_exceedanceProbability: cannot compute analytical ep for non diagonal covariance. Please use the sampling method');
                end
                dSigma = max (diag (Sigma), 1e-6);
                mu = VBA_vec(mu);
                for k = 1 : K
                    l = setdiff (1 : K, k);  
                    % p(x_k > x_j) = int pdf_k(t)cdf_1(t)cdf_2(t)...cdf_k-1(t)cdf_k+1(t)...cdf_k(t) dt 
                    mycdf = @(t) sum (bsxfun(@(x,i) log(VBA_spm_Ncdf(x, mu(l(i)), dSigma(l(i)))), VBA_vec(t), 1 : (K-1)), 2)';                 
                    ii = @(t) exp(log(VBA_spm_Npdf(t, mu(k), dSigma(k))) + mycdf(t));
                    ep(k) = integral(ii, -Inf, Inf);
                end
                ep = ep / sum(ep);
        end
    % ---------------------------------------------------------------------
    case 'Dirichlet'
        
        switch options.method
            
            case 'sampling'
                r_samp = VBA_random ('Dirichlet', alpha, options.nSamples);
                
                ep = nansum (bsxfun (@eq, r_samp, max (r_samp)), 2);
                ep = ep / nansum (ep);
                
                if VBA_isWeird (r_samp) 
                    warning ('VBA_exceedanceProbability: unstable parametrization (alpha < 1), think about trying the analytical form.');
                end
                 
            case 'analytical'
                alpha = min (alpha, 1e7);
                
                for k = 1 : K
                    f = @(x) integrand (x, alpha(k), alpha(1 : K ~= k));
                    ep(k) = integral (f, eps, Inf, 'Waypoints', 10.^(-15:15));
                end
                ep = ep / sum (ep);
        end
     
end

% formating
% =========================================================================
ep = VBA_vec (ep);

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

