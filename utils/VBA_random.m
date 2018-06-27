function X = VBA_random (name, varargin)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% X = VBA_random (name, p1, p2, ..., [N, M, ...])
% generate random numbers following the requested distribution
%
% If you the MATLAB statistics toolbox, this merely serve as an overload of 
% the random() function, although with a slightly different syntax and
% behaviour. If the statistics toolbox is not installed, VBA_random will 
% fall back to in-house generators (mainly using SPM routines).
% 
%
% IN:
%   - name: type of distribution (see below)
%   - p1, p2, ...: parameters of the distribution
%   - N, M, ...: dimensions of the output. Default = 1 sample.
%       + for univariate densities, returns a N x N array if only one
%         dimension is provided, or a N x M x ... array if more dimensions 
%         are specified.
%       + for multivariate densities, samples are returned as columns of X 
%         and only N can be used to specify the number of samples requested
%
% OUT:
%   - X: random samples
%
% Available distributions:
% ~~~~~~~~~~~~~~~~~~~~~~~~
%
% 'Arbirary'
% ----------
%     + parameters:
%         - p: k x 1 vector  describing the density of possible values
%         - vals: k x 1 vector or k x j matrix of values
%     + output:
%         - each element (if val is a vector) or column (if val is an array)
%           of X will follow the distribution P(x = val(v)) = p(v), or 
%           P(x = val(v, :)) = p(v) respectively.
%
% 'Bernoulli'
% -----------
%     + parameter:
%        - p: parameter of the Bernoulli distribution
%     + output: samples from the Bernoulli distribution, ie such that 
%           P(x = 1) = p and P(x = 0) = 1 - p
%
% 'Binomial'
% ----------
%     + parameters:
%           - n: number of trials
%           - p: probability of success of each trial
%     + output: samples from the distribution B(n, p), defined as the sum 
%           of n independent samples from a Bernoulli distribution with
%           a parameter p
%
% 'Categorical'
% -------------
%     + parameters:
%           - p: vector defining the probability of the categories
%     + output: integer k such that P(x = k) = p(k). Note that Categorical
%           samples are sometimes expressed as vectors with all elements
%           set to 0 except for the winning category element set to 1. 
%           Use VBA_indicator (X, k) for the conversion.
%
% 'Dirichlet'
% -----------
%     + parameters:
%           - alphas: pseudo-counts
%     + output: points of the simplex (ie vectors summing to 1) sampled
%           from the Dirichlet distribution D(alphas)
%
% 'Gamma'
% -------
%     + parameters:
%           - a: shape parameter
%           - b: scale parameter
%     + output: samples from the Gamma distribution Ga(a, b). Note that the
%         corresponding expectation is therefore a * b (not a / b).
%
% 'Gaussian'
% ----------
%     + parameters:
%           - mu: scalar or k x 1 vector, mean of the Normal distribution
%           - Sigma2: scalar or k x k array, (co-)variance of
%     + output: if parameters are scalars, samples from the univariate
%           normal distribution N(mu, Sigma).
%           If mu and Sigma are respectively a vector and an array of
%           matching dimensions, columns of X are samples from the 
%           corresponding multivariate Gaussian distribution.
%
% 'Multinomial'
% -------------
%     + parameters:
%           - n: number of trials
%           - p: 1 x k vector, probability of each category winning a trial
%     + output: each column of X is a k-vector that sums to n. It is
%           constructed as the sum of n Categorical samples expressed in 
%           their indicator vector form.
%
% /////////////////////////////////////////////////////////////////////////

    switch name
        case 'Arbitrary'
            % get parameters
            [p, val, N] = getParam (varargin);
            assert (isvector (p) && numel (p) > 1 && abs(sum (p) - 1) < 1e-15, ...
                'VBA:invalidInput', ...
                '*** VBA_random: p must be a vector summing to one.');
            if isvector (val)
                assert (numel (val) == numel (p), ...
                    'VBA:invalidInput', ...
                    '*** VBA_random: inconsistent sizes');
                val = VBA_vec (val);
                if isscalar (N)
                    N = {N{1}, N{1}};
                end
            else
                assert (size (val, 1) == numel (p), ...
                    'VBA:invalidInput', ...
                    '*** VBA_random: inconsistent sizes');
                assert (isscalar (N), ...
                    'VBA:invalidInput', ...
                    '*** VBA_random: N must be scalar for multivariate sampling.');
            end
            % sample
            pcdfu = cumsum (p(:));
            pcdfl = [0; pcdfu(1 : end - 1)];
            
            s = rand (1, prod ([N{:}]));
            [idx, ~] = find (bsxfun (@le, s, pcdfu) & bsxfun (@gt, s, pcdfl));
            
            X = val(idx,:)';
            
            if isvector (val)     
                X = reshape (X, N{:});
            else
                assert (isscalar(N), '*** VBA_random: N must be scalar for multivariate sampling.');
            end
            
        case 'Bernoulli'
            % get parameters
            [p, N] = getParam (varargin);
            assert (VBA_isInRange (p(~ isnan(p)), [0, 1]), ...
                'VBA:invalidInput', ...
                '*** VBA_random: p must be between 0 and 1.');
            if ~ isscalar (p) 
                if isscalar (N) && N{1} == 1
                    N = num2cell (size (p));
                else
                    assert (all (size (p) == [N{:}]), ...
                    'VBA:invalidInput', ...
                    '*** VBA_random: size information is inconsistent.');
                end
            else
                if isscalar (N)
                    N = {N{1}, N{1}};

                end
            end
            % sample
            try
                X = binornd (1, p, N{:});
            catch
                X = + (rand (N{:}) <= p);
            end
            X(isnan (p)) = nan;
            
        case 'Binomial'
            % get parameters
            [n, p, N] = getParam (varargin);
            assert (isscalar (n) && n > 0 && rem (n, 1) == 0, ...
                'VBA:invalidInput', ...
                '*** VBA_random: n must be a positive integer.');
            assert (isscalar (p) && VBA_isInRange (p, [0, 1]), ...
                'VBA:invalidInput', ...
                '*** VBA_random: p must be a scalar between 0 and 1.');
            % sample
            try
                X = binornd (n, p, N{:});
            catch 
                N{end + 1} = n;
                X = sum (VBA_random ('Bernoulli', p, N{:}), numel(N));
            end 
            
        case 'Categorical'
            [p, N] = getParam (varargin); 
            X = VBA_random ('Arbitrary', p, 1 : numel(p), N{:});
      
        case 'Dirichlet'
            % get parameters
            [alpha, N] = getParam (varargin); 
            assert (isvector (alpha) && numel (alpha) > 1 && all (alpha > 0), ...
                'VBA:invalidInput', ...
                '*** VBA_random: alpha must be a vector of positive values.');
            assert (isscalar (N), ...
                'VBA:invalidInput', ...
                '*** VBA_random: N must be scalar for multivariate sampling.');
            N = N{1};
            % sample
            K = numel (alpha);
            scale = 1;
            try
                r = gamrnd (repmat(alpha(:), 1, N), scale, K, N);
                r(isinf (r)) = realmax;
            catch 
                r = zeros (K, N);
                for k = 1 : K
                    r(k, :) = VBA_spm_gamrnd (alpha(k), scale, N);
                end
            end
            X = bsxfun (@rdivide, r, sum(r));
           
        case 'Gamma'
             % get parameters
            [a, b, N] = getParam (varargin);  
            assert (isscalar (a) && isscalar(b) && all ([a, b] > 0), ...
                'VBA:invalidInput', ...
                '*** VBA_random: a and b must be positive scalars.');
            % sample
            try
                X = gamrnd (a, b, N{:});
            catch
                X = VBA_spm_gamrnd (a, b, N{:});
            end
            
        case 'Gaussian'
            [mu, Sigma, N] = getParam (varargin);
            assert (all (size (Sigma) == numel (mu)), ...
                'VBA:invalidInput', ...
                '*** VBA_random: inconsistent sizes');

            % + univariate
            if isscalar (mu)
                if isscalar (N)
                    N = {N{1}, N{1}};
                end
                X = mu + sqrt (Sigma) * randn (N{:}); 
                
            % + multivariate
            else
                assert (isscalar (N), ...
                    'VBA:invalidInput', ...
                    '*** VBA_random: N must be scalar for multivariate sampling.');
                assert (VBA_issymmetric (Sigma), ...
                    'VBA:invalidInput', ...
                    '*** VBA_random: Sigma must be symmetric positive definite.');
                mu = VBA_vec (mu);
                k = numel (mu);
                X = bsxfun (@plus, VBA_sqrtm (Sigma) * randn (k, N{1}), mu) ;  
            end        
            
            
        case 'Multinomial'
            % get parameters
            [n, p, N] = getParam (varargin); 

            assert (isscalar (N), ...
                 'VBA:invalidInput', ...
                 '*** VBA_random: N must be scalar for multivariate sampling.');
            if any (isnan (p))
                X = nan (numel(p), N{1});
                return
            end
            assert (isscalar (n) && n > 0 && rem (n, 1) == 0, ...
                'VBA:invalidInput', ...
                '*** VBA_random: n must be a positive integer.');
            assert (isvector (p) && numel (p) > 1 && abs(sum (p) - 1) < 1e-15, ...
                'VBA:invalidInput', ...
                '*** VBA_random: p must be a vector summing to one.');
            % sample
            X = zeros (numel (p), N{1});
            for i = 1 : n
                X = X + VBA_random ('Arbitrary', p, eye (numel(p)), N{1});
            end
            
        otherwise
            error('*** VBA_random: Unkown distribution type');
    end
end

function varargout = getParam (args)
    varargout = args(1 : nargout - 1);
    dims = args(nargout : end);
    if isempty (dims)
        dims = {1};
    end
    varargout{nargout} = dims;
end