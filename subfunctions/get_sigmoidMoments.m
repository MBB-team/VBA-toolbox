function [ms,vs] = get_sigmoidMoments(muX,VX)

% Laplace approximation of Gaussian processes passed through the sigmoid
% function [ms,vs] = get_sigmoidMoments(muX,VX)
% IN:
%   - muX: nxt time series of 1st-order moments
%   - VX: nxt time series of 2nd-order moments
% OUT:
%   - ms,vs: nxt time series of 1st- and 2nd-order moments of the mapped
%   process


g_fname = @g_sigm;
opt.inG.x = 0;
opt.priors.a_sigma = Inf;
opt.priors.b_sigma = 0;

dim.n = 0;
dim.n_theta = 0;
dim.n_phi = size(muX,1);
dim.n_t = 1;
dim.p = size(muX,1);

n_t = size(muX,2);
ms = zeros(size(muX));
vs = zeros(size(VX));
for t=1:n_t
    opt.priors.muPhi = muX(:,t);
    opt.priors.SigmaPhi = diag(VX(:,t));
    [ms(:,t),VS] = VBA_getLaplace(0,[],g_fname,dim,opt);
    vs(:,t) = diag(VS);
end



