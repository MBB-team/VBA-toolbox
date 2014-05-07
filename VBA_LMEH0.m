function [LLH0] = VBA_LMEH0(y,options)
% evaluates the log-evidence of the null (H0) model
% [LLH0] = VBA_LMEH0(y)
% This function evaluates the exact log marginal likelihood of a null model
% H0, which basically assumes that the data only contain Gaussian noise of
% unknown variance, i.e. y ~ N(0,sigma^-1*I), where sigma is unknown, with
% Jeffrey's priors, i.e. sigma ~ Gamma(1,1). This is a trivial (null)
% model, to be compared with the model evidence of a non-trivial model as a
% diagnostic measure.
% IN:
%   - y: the data matrix;
% OUT:
%   - LLH0: the log evidence of the null model

if isfield(options,'extended') && options.extended
    LLH0=VBA_LMEH0_extended(y,options);
    return;
end

try
    y = y(options.isYout==0);
end
n = numel(y);
if ~isbinary(y)
    try
        a0 = options.priors.a_sigma;
        b0 = options.priors.b_sigma;
    catch
        a0 = 1;
        b0 = 1;
    end
    y = y(:);
    y2 = y'*y;
    alpha = n/2 + a0;
    beta = y2./2 + b0;
    LLH0 = -0.5*n*log(2*pi) + a0.*log(b0) - alpha.*log(beta) - gammaln(a0) + gammaln(alpha);
else % assume binary data under chance (p=0.5)
    LLH0 = n*log(.5);
end
    

