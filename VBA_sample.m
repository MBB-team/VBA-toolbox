function y = VBA_sample(form,suffStat,N,verbose)
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

if verbose
    fprintf(1,['Sampling from ',form,' distribution... ']);
    fprintf(1,'%6.2f %%',0)
end
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
                if mod(i,N./20) < 1 && verbose
                    fprintf(1,repmat('\b',1,8))
                    fprintf(1,'%6.2f %%',100*i/N)
                end
            end
        end
        
    case 'dirichlet'
        K = size(suffStat.d,1);
        try
            r = gamrnd(repmat(vec(suffStat.d),1,N),1,K,N);
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
                if mod(i,N./20) < 1 && verbose
                    fprintf(1,repmat('\b',1,8))
                    fprintf(1,'%6.2f %%',100*i/N)
                end
            end
        end
        
    case 'multinomial'
        try
            y = mnrnd(suffStat.n,suffStat.p,N)';
        catch
            K = size(suffStat.p,1);
            y = zeros(K,N);
            for i=1:suffStat.n
                y = y + sampleFromArbitraryP(suffStat.p,eye(K),N)';
                if verbose
                    fprintf(1,repmat('\b',1,8))
                    fprintf(1,'%6.2f %%',100*i/suffStat.n)
                end
            end
        end
        
end
if verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
end

function r = drchrnd(a,n)
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

