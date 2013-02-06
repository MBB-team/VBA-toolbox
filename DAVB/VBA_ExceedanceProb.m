function ep = VBA_ExceedanceProb(mu,Sigma,form,verbose)
% calculates the exceedance probability for mutivariate Gaussian variables
% function ep = VBA_ExceedanceProb(mu,Sigma)
% IN:
%   - mu/Sigma: first- and second-order moments of the Gaussian density
%   - form: 'gaussian' or 'dirichlet'
% OUT:
%   - ep: vector of exceedance probabilities, i.e. the probability, for
%   each variable, to be greater than all the other ones.

try, form; catch form='gaussian'; end
try, verbose; catch verbose=0; end

K = size(mu,1);
ep = ones(K,1);
c = [1;-1];
switch form
    case 'gaussian'
        for k=1:K
            for l=setdiff(1:K,k)
                ind = [k,l];
                m = mu(ind);
                V = Sigma(ind,ind);
                ep(k) = ep(k)*VB_PPM(c'*m,c'*V*c,0,0);
            end
        end
        ep = ep./sum(ep);
    case 'dirichlet'
        ep = myget_ep(mu,1e4,verbose);
end



function xp = myget_ep(alpha,Nsamp,verbose)
% sampling approx for dirichlet disitribution
if verbose
    fprintf(1,'Evaluating exceedance probability...');
    fprintf(1,'%6.2f %%',0)
end
Nk = size(alpha,1);
r_samp = zeros(Nsamp,Nk);
for samp=1:Nsamp
    for k = 1:Nk
        r(:,k) = spm_gamrnd(alpha(k),1);
    end
    sr = sum(r,2);
    for k = 1:Nk
        r(:,k) = r(:,k)./sr;
    end
    r_samp(samp,:)=r;
    if mod(samp,1e3) < 1 && verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',100*samp/Nsamp)
    end
end
xp = zeros(1,Nk);
[y,j]=max(r_samp,[],2);
tmp=histc(j,1:Nk)';
xp=tmp/Nsamp;
if verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
end
