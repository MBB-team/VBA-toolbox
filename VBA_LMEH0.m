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


LLH0 = 0;

% gaussian source
s_g = find([options.sources(:).type]==0);
for s=1:length(s_g)
    y_i = options.sources(s_g(s)).out ;
    y_s = y(y_i,:);
    y_s = y_s(options.isYout(y_i,:)==0);
    LLH0 = LLH0 + lev_GLM(VBA_vec(y_s),ones(numel(y_s),1));
%     try
%         a0 = options.priors.a_sigma(s);
%         b0 = options.priors.b_sigma(s);
%     catch
%         a0 = 1;
%         b0 = 1;
%     end
%     n = numel(y_s);
%     y_s = y_s(:);
%     y2 = y_s'*y_s;
%     alpha = n/2 + a0;
%     beta = y2./2 + b0;
%     LLH0 = LLH0 -0.5*n*log(2*pi) + a0.*log(b0) - alpha.*log(beta) - gammaln(a0) + gammaln(alpha);
end

% binomial sources
s_b = find([options.sources(:).type]==1);
for s=1:length(s_b)
    y_i = options.sources(s_b(s)).out ;
    y_s = y(y_i,:);
    y_s = y_s(options.isYout(y_i,:)==0);
    n = numel(y_s);
    LLH0 = LLH0 + n*log(.5);
end

% multinomial sources
s_m = find([options.sources(:).type]==2);
for s=1:length(s_m)
    y_i = options.sources(s_m(s)).out ;
    y_s = y(y_i,:);
    isYinIdx = find(~any(options.isYout(y_i,:)));
    y_s = y_s(:,isYinIdx);
    [c,n] = size(y_s); 
    LLH0 = LLH0 + n*log(1/c);
end
