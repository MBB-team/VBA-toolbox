function ep = VBA_ExceedanceProb(mu,Sigma)
% calculates the exceedance probability for mutivariate Gaussian variables
% function ep = VBA_ExceedanceProb(mu,Sigma)
% IN:
%   - mu/Sigma: first- and second-order moments of the Gaussian density
% OUT:
%   - ep: vector of exceedance probabilities, i.e. the probability, for
%   each variable, to be greater than all the other ones.
K = size(mu,1);
ep = ones(K,1);
c = [1;-1];
for k=1:K
    for l=setdiff(1:K,k)
        ind = [k,l];
        m = mu(ind);
        V = Sigma(ind,ind);
        ep(k) = ep(k)*VB_PPM(c'*m,c'*V*c,0,0);
    end
end
ep = ep./sum(ep);