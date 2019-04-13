function H = VBA_entropy (name, varargin)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% H = VBA_entropy (name, p1, p2, ...)
% compute the entropy of the requested distribution 
%
% IN:
%   - name: type of distribution (see below)
%   - p1, p2, ...: parameters of the distribution
%
% OUT:
%   - H: entropy
%
% Available distributions:
% ~~~~~~~~~~~~~~~~~~~~~~~~
%
% 'Bernoulli'
% -----------
%     + parameter:
%        - p: parameter of the Bernoulli distribution
%     + output: 
%
%           H = - q log(q) - (1 - p) log(1 - p)
%
% 'Binomial'
% ----------
%     + parameters:
%           - n: number of trials
%           - p: probability of success of each trial
%     + output: 
%                          |n|                       |n|
%           H = sum(k=0,n) | | p^k (1-p)^(n-k) log [ | | p^k (1-p)^(n-k) ] 
%                          |k|                       |k|
% 'Categorical'
% -------------
%     + parameters:
%           - p: vector defining the probability of the categories
%     + output: 
%
%           H = - sum p_i log(p_i) 
%
% 'Dirichlet'
% -----------
%     + parameters:
%           - alphas: vector of K pseudo-counts
%     + output: 
%           with alpha_0 = sum(alphas)
%
%           H = log Beta(alphas) + (alphas0 - K) psi(alpha0) - sum_j (alpha_j - 1) psi(alpha_j) 
%
% 'Gamma'
% -------
%     + parameters:
%           - theta: shape parameter
%           - k: scale parameter 
%     + output: 
%
%           H = k + log theta + log Gamma(k) + (1 - k) psi(k)
%
%     Beware the notation. k is the scale parameter, i.e. the expecation is
%     therefore theta * k (not theta / k). For the rate notation, use
%     V_entropy ('Gamma', shape, 1/rate).
%
% 'Gaussian'
% ----------
%     + parameters:
%           - mu: scalar or 1 x k vector, mean of the distribution
%           (optional)
%           - Sigma2: scalar or k x k array, (co-)variance of the  distribution
%     + output:
%
%           H = 1/2 log det (2 pi e Sigma2)
%
%     As the entropy only depends on the (co-)variance, you can call the 
%     shortcut VBA_entropy ('Gaussian', Sigma2) instead of 
%     VBA_entropy ('Gaussian', mu,  Sigma2).
%
% 'Multinomial'
% -------------
%     + parameters:
%           - n: number of trials
%           - p: 1 x k vector, probability of each category winning a trial
%     + output: 
%
%           H = - log(n!) - n sum [p_i log(p_i)]
%                                         | n |
%               + sum(i=1,k) sum(x_i=0,n) |   | p_i^x_i (1 - p_i)^(n - x_i) log(x_i!)
%                                         |x_i|
%
% /////////////////////////////////////////////////////////////////////////

switch name
    
    case 'Bernoulli'
        % get parameters
        p = varargin{1};
        % entropy
        q = 1 - p;
        H = - q * log (q) - p * log (p);
        
    case 'Binomial'
        % get parameters
        [n, p] = varargin{:};
        % entropy
        H = 0;
        for k = 0 : n
            H_i = nchoosek(n,k) * p^k * (1 - p)^(n - k);
            H = H - H_i * log (H_i);
        end
        
    case 'Categorical'
        % get parameters
        p = [varargin{:}];
        % entropy
        H = - sum (p .* log (p));
        
    case 'Dirichlet'
        % get parametesr
        alphas = [varargin{:}];
        % entropy
        K = numel (alphas);
        alpha_0 = sum (alphas);
        log_B_alpha = sum (gammaln (alphas)) - gammaln (sum (alphas));
        H = log_B_alpha + (alpha_0 - K) * psi (alpha_0) - sum ((alphas - 1) .* psi (alphas));
        
    case 'Gamma'
        % get parameters 
        [theta, k] = varargin{:};
        % entropy
        H = theta + log (k) + gammaln (theta) + (1 - theta) .* psi (theta);
        
    case 'Gaussian'
        % get parameters 
        Sigma2 = varargin{end};
        % entropy
        n = size (Sigma2, 1);
        H = 0.5 * (n * (1 + log (2 * pi)) + VBA_logDet (Sigma2));
        
    case 'Multinomial'
        % get parameters
        [n, p] = varargin{:};
        % entropy
        H = - log (factorial (n)) - n * sum (p .* log (p));
        for xi = 0 : n
            H = H + sum (nchoosek (n, xi) * p.^xi .* (1 - p).^(n-xi) * log (factorial (xi)));
        end

    otherwise
        error ('*** VBA_entropy: unkown distribution type');
end

