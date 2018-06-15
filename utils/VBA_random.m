function X = VBA_random (name, varargin)

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
            assert (isscalar (p) && VBA_isInRange (p, [0, 1]), ...
                'VBA:invalidInput', ...
                '*** VBA_random: p must be a scalar between 0 and 1.');
            % sample
            try
                X = binornd (1, p, N{:});
            catch
                X = + (rand (N{:}) <= p);
            end
            
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
                mu = VBA_vec (mu);
                k = numel (mu);
                X = bsxfun (@plus, sqrtm (Sigma) * randn (k, N{1}), mu) ;  
            end        
            
            
        case 'Multinomial'
            % get parameters
            [n, p, N] = getParam (varargin); 
            assert (isscalar (n) && n > 0 && rem (n, 1) == 0, ...
                'VBA:invalidInput', ...
                '*** VBA_random: n must be a positive integer.');
            assert (isvector (p) && numel (p) > 1 && abs(sum (p) - 1) < 1e-15, ...
                'VBA:invalidInput', ...
                '*** VBA_random: p must be a vector summing to one.');
            assert (isscalar (N), ...
                 'VBA:invalidInput', ...
                 '*** VBA_random: N must be scalar for multivariate sampling.');
                
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
    N = args(nargout:end);
    if isempty (N)
        N = {1};
    end
    varargout{nargout} = N;
end