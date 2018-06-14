function X = VBA_random (name, varargin)


    switch name
        case 'Arbitrary'
            
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
            
        case 'Dirichlet'
            
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
            
        case 'Multinomial'
            
            
        otherwise
    end
end

function varargout = getParam (args)
    for i = 1 : nargout - 1
        varargout{i} = args{i};
    end
    N = args(nargout:end);
    if isempty (N)
        N = {1};
    end
    varargout{nargout} = N;
end