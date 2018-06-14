function y = VBA_sample(form,suffStat,N)
% samples from exponential family probability distribution functions
% function y = VBA_sample(form,ss,N,verbose)
% IN:
%   - form: 'gaussian', 'gamma', 'dirichlet', 'multinomial' or '1D'
%   - suffStat: a structure with appropriate sufficient statistics, i.e.:
%       -> if form='gaussian', suffStat.mu = E[y] and suffStat.Sigma = V[y]
%       -> if form='gamma', suffStat.a = shape parameter, and suffStat.b =
%       scale parameter of the Gamma density
%       -> if form='dirichlet', suffStat.d = Dirichlet counts
%       -> if form='multinomial', suffStat.p = multinomial probabilities
%       and suffstat.n = number of independent trials
%   - N: number of samples
%   - verbose: verbose mode
% OUT:
%   - y: KXN array of vector-valued samples (where K is the dimension of
%   the sampled data).
% NOTE: by default, this function tries to use Matlab pseudo-random
% samplers. It reverts to SPM in case these functions cannot be called.

switch form
    
    case 'gaussian'
        S = VBA_getISqrtMat(suffStat.Sigma,0);
        
        n = size(suffStat.mu,1);
        y = repmat(suffStat.mu,1,N) + S*randn(n,N);
        
    case 'gamma'
        try
            y=gamrnd(suffStat.a,suffStat.b,1,N);
        catch
            y = zeros(1,N);
            for i=1:N
                y(i) = VBA_spm_gamrnd(suffStat.a,suffStat.b);
            end
        end
        
    case 'dirichlet'
        K = size(suffStat.d,1);
        try
            r = gamrnd(repmat(VBA_vec(suffStat.d),1,N),1,K,N);
            r(isinf(r)) = realmax;
            y = r ./ repmat(sum(r,1),K,1);
        catch
            y = zeros(K,N);
            r = zeros(K,1);
            for i=1:N
                for k = 1:K
                    r(k) = VBA_spm_gamrnd(suffStat.d(k),1);
                end
                y(:,i) = r./sum(r);
            end
        end
        
    case 'bernoulli'
        try
            y = binornd (1, suffStat.p, 1, N);
        catch 
            y = + (rand (1, N) <= suffStat.p);
        end
        
    case 'binomial'
        try
            y = binornd (suffStat.n, suffStat.p, 1, N);
        catch 
            y = zeros (1, N);
            for i = 1 : suffStat.n
                y = y + VBA_sample ('bernoulli', suffStat, N);
            end
        end
        
    case 'multinomial'
        try
            y = mnrnd(suffStat.n,suffStat.p,N)';
        catch
            K = numel(suffStat.p);
            y = zeros(K,N);
            for i=1:suffStat.n
                y = y + VBA_indicator(VBA_sampleFromArbitraryP(suffStat.p, 1:K, N), K);
            end
        end
        
end
