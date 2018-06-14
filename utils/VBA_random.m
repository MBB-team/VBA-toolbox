function X = VBA_random (name, varargin)


    switch name
        case 'Arbitrary'
            % get parameters
            [p, val, N] = getParam (varargin);
            if isvector (val)
                val = VBA_vec (val);
                if isscalar (N)
                    N = {N{1}, N{1}};
                end
            end
            % sample
            pcdfu = cumsum (p(:));
            pcdfl = [0; pcdfu(1 : end - 1)];
            
            s = rand (1, prod ([N{:}]));
            [idx, ~] = find (s <= pcdfu & s > pcdfl);
            
            X = val(idx,:)';
            
            if isvector (val)     
                X = reshape (X, N{:});
            else
                assert (isscalar(N), '*** VBA_random: N must be scalar for multivariate sampling.');
            end
            
        case 'Bernoulli'
            % get parameters
            [p, N] = getParam (varargin);
            % sample
            try
                X = binornd (1, p, N{:});
            catch
                X = + (rand (N{:}) <= p);
            end
            
        case 'Binomial'
            % get parameters
            [n, p, N] = getParam (varargin);       
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
            [alphas, N] = getParam (varargin); 
            assert (numel (N) == 1, '*** VBA_random: N must be scalar for multivariate sampling.');
            N = N{1};
            % sample
            K = numel (alphas);
            scale = 1;
            try
                r = gamrnd (repmat(alphas(:), 1, N), scale, K, N);
                r(isinf (r)) = realmax;
            catch 
                r = zeros (K, N);
                for k = 1 : K
                    r(k, :) = VBA_spm_gamrnd (alphas(k), scale, N);
                end
            end
            X = bsxfun (@rdivide, r, sum(r));
           
        case 'Gamma'
             % get parameters
            [a, b, N] = getParam (varargin);       
            % sample
            try
                X = gamrnd (a, b, N{:});
            catch
                X = VBA_spm_gamrnd (a, b, N{:});
            end
            
        case 'Gaussian'
            [mu, Sigma, N] = getParam (varargin);
            if ~ isscalar(mu) || ~ isscalar(Sigma)
                assert(isscalar (N), '*** VBA_random: N must be scalar for multivariate sampling.');
                mu = VBA_vec(mu);
                assert(all(size(Sigma) == numel(mu)),'*** VBA_random: inconsistent moments size');
            else
                if isscalar (N)
                    N = {N{1}, N{1}};
                end
            end        
            n = numel(mu);
            X = bsxfun(@plus, randn (N{:}, n) * sqrtm(Sigma) , mu')' ;
            
        case 'Multinomial'
            % get parameters
            [n, p, N] = getParam (varargin); 
            assert(isscalar (N), '*** VBA_random: N must be scalar for multivariate sampling.');
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