function ep = VBA_ExceedanceProb(mu,Sigma,form,verbose,Nsamp)
% calculates the exceedance probability for mutivariate Gaussian variables
% function ep = VBA_ExceedanceProb(mu,Sigma)
% IN:
%   - mu/Sigma: sufficient statistics of the pdf
%       -> if form='gaussian': mu=E[x] and Sigma=V[x]
%       -> if form='dirichlet': mu=Dirichlet counts and Sigma is unused
%   - form: 'gaussian' or 'dirichlet'
% OUT:
%   - ep: vector of exceedance probabilities, i.e. the probability, for
%   each variable, to be greater than all the other ones.

try, form; catch, form = 'gaussian'; end 
try, verbose; catch, verbose=0; end
try, Nsamp; catch, Nsamp=1e4; end

K = size(mu,1);
ep = ones(K,1);
c = [1;-1];
switch form
    case 'gaussian'
        r_samp = VBA_sample('gaussian',struct('mu',mu,'Sigma',Sigma),Nsamp,verbose)';
        [y,j] = max(r_samp,[],2);
        tmp = histc(j,1:length(mu))';
        ep = tmp/Nsamp;
    case 'gaussian2'
        for k=1:K
            for l=setdiff(1:K,k)
                ind = [k,l];
                m = mu(ind);
                V = Sigma(ind,ind);
                ep(k) = ep(k)*VBA_PPM(c'*m,c'*V*c,0,0);
            end
        end
        ep = ep./sum(ep);
    case 'dirichlet'
        r_samp = VBA_sample('dirichlet',struct('d',mu),Nsamp,verbose)';
        if isweird(r_samp)
            error('*** Could not compute the exceedance probability because of numerical instabilities.');
        end
        [y,j] = max(r_samp,[],2);
        tmp = histc(j,1:length(mu))';
        ep = tmp/Nsamp;
end

