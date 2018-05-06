function [fx] = f_VBvolatile0(x,P,u,in)
% computes Laplace-VB update rule of HGF learner
% [fx] = f_VBvolatile0(x,P,u,in)
% This is the one-step Markovian update rule for the posterior sufficient
% statistics of a volatile environment, as derived in [Mathys et al. 2010].
% Note: HGF stands for "Hierarchical Gaussian Filter".
% IN:
%   - x: the previous posterior sufficient statistics:
%   x(1)= u [U is the outcome, wose probability is tracked]
%   x(2)= E[log-odds of P(u=1)]
%   x(3)= log V[log-odds of P(u=1)]
%   x(4)= E[log-volatility]
%   x(5)= log V[log-volatility]
%   - P: the perceptual model parameters vector, ie. P = [ka;om;th], using
%   the notation of [Mathys et al. 2010].
%   - u: the outcome, whose probability is tracked over trials.
%   - in: options set in options.inF
% OUT:
%   - fx: the updated posterior sufficient statistics (having observed u).

% transform and define states and parameters
x(3) = exp(x(3)); % variance on second-level states is in log-space
x(5) = exp(x(5)); % variance on third-level states is in log-space
ka = in.lev2 * VBA_sigmoid(P(1), 'scale', in.kaub);
om = P(2); % volatility rescaling
th = VBA_sigmoid(P(3), 'scale', in.thub);
vol = exp(ka*x(4)+om);
fx = zeros(size(x));

% 1st level (perceptual level --> degrade if subliminal!)
fx(1) = u(1); % trivial first-level states

% 2nd level
s1h = VBA_sigmoid(x(2))*(1-VBA_sigmoid(x(2))); % likelihood precision
pe1 = fx(1) - VBA_sigmoid(x(2)); % prediction error
s2h = x(3) + vol; % 2nd-level prediction variance
fx(3) = 1/(s2h^-1 + s1h); % posterior variance
fx(2) = x(2) + fx(3)*pe1; % 2nd-level update
fx(2) = VBA_sigmoid(VBA_sigmoid(fx(2)),'inverse',true); % for numerical purposes

% 3rd level
pi3h = 1/(x(5)+th);
w2 = 1./(1+x(3)/vol);
r2 = (1-x(3)/vol)./(1+x(3)/vol);
pe2 = (fx(3)+(fx(2)-x(2))^2)/(vol+x(3)) -1;
pl = max([0,.5*ka^2*w2*(w2+r2*pe2)]); % for numerical purposes
fx(5) = 1/(pi3h + pl);
fx(4) = x(4) + .5*ka*w2*fx(5)*pe2;

% retransform states
fx(3) = log(fx(3));
fx(5) = log(fx(5));

