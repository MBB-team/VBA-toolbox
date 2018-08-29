function y = VBA_sample(form,suffStat,N)
% legacy code
s = warning ('on');
warning ('*** The function `VBA_sample` is now deprecated. Please see `VBA_random` for an alternative.') 
warning (s);

if nargin < 3
    N = 1;
end

switch form
    case 'gaussian'
        if isscalar(suffstat.mu)
            N = {1, N};
        else
            N = {N};
        end
        
        y = VBA_random ('Gaussian', suffStat.mu, suffStat.Sigma, N{:});
        
    case 'gamma'
        y = VBA_random ('Gamma', suffStat.a, suffStat.b, 1, N);
    
    case 'dirichlet'
        y = VBA_random ('Dirichlet', suffStat.d, N);
        
     case 'bernoulli'
        y = VBA_random ('Bernoulli', suffStat.p, 1, N);
        
    case 'binomial'
        y = VBA_random ('Binomial', suffStat.n, suffstat.p, 1, N);
        
    case 'multinomial'
        y = VBA_random ('Multinomial', suffStat.n, suffstat.p, N);

end
