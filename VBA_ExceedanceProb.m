function ep = VBA_ExceedanceProb(mu,Sigma,form,verbose,Nsamp)
% calculates the exceedance probability for mutivariate variables
% function ep = VBA_ExceedanceProb(mu,Sigm,forma)
% IN:
%   - mu/Sigma: sufficient statistics of the pdf
%       -> if form = 'gaussian' : mu = E[x] and Sigma = V[x]
%       -> if form = 'dirichlet': mu = Dirichlet alphas and Sigma is unused
%   - form: 'gaussian', 'gaussian2', 'dirichlet' or 'dirichlet2'
%       -> if form = 'gaussian'  : sampling solution for mvn
%       -> if form = 'gaussian2' : analytical solution for mvn
%       -> if form = 'dirichlet' : sampling solution for Dirichlet
%       -> if form = 'dirichlet2': analytical solution for Dirichlet
%   - verbose: verbose mode for random sampling
%   - Nsamp: number of samples for random sampling
% OUT:
%   - ep: vector of exceedance probabilities, i.e. the probability, for
%         each variable, to be greater than all the other ones.
% NOTE: The analytical solution for the Dirichlet distribution is based
% on numerical integration over Gamma distributions [1,2] and was added
% by Joram Soch [3] on 17/05/2018.
% [1] https://arxiv.org/abs/1611.01439
% [2] https://github.com/JoramSoch/MACS/blob/master/MD_Dir_exc_prob.m
% [3] mailto:joram.soch@bccn-berlin.de

try, form; catch, form = 'gaussian'; end 
try, verbose; catch, verbose = 0; end
try, verbose; catch, Nsamp = 1e4; end

K  = size(mu,1);
ep = ones(K,1);

switch form
    
    case 'gaussian'
        r_samp = VBA_sample('gaussian',struct('mu',mu,'Sigma',Sigma),Nsamp,verbose)';
        [y,j]  = max(r_samp,[],2);
        tmp    = histc(j,1:length(mu))';
        ep     = tmp/Nsamp;
    
    case 'gaussian2'
        c = [1; -1];
        for k = 1:K
            for l = setdiff(1:K,k)
                ind = [k,l];
                m   = mu(ind);
                V   = Sigma(ind,ind);
                ep(k) = ep(k)*VBA_PPM(c'*m,c'*V*c,0,'gaussian',0);
            end
        end
        ep = ep./sum(ep);
    
    case 'dirichlet'
        r_samp = VBA_sample('dirichlet',struct('d',mu),Nsamp,verbose)';
        [y,j]  = max(r_samp,[],2);
        tmp    = histc(j,1:length(mu))';
        ep     = tmp/Nsamp;
    
    case 'dirichlet2'
        alpha = mu;
        for k = 1:K
            f = @(x) integrand(x,alpha(k),alpha([1:K]~=k));
            ep(k) = integral(f,eps,Inf);
        end
        ep = ep./sum(ep);
        
end

% Integrand function for numerical integration
%-------------------------------------------------------------------------%
function p = integrand(x,aj,ak)

p = ones(size(x));
for k = 1:numel(ak)
    p = p .* gammainc(x,ak(k));                         % Gamma CDF
end
p = p .* exp((aj-1).*log(x) - x - gammaln(aj));         % Gamma PDF