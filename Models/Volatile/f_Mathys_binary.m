function [fx] = f_Mathys_binary(x,P,u,in)
% Compute the VB-laplace update rules for hidden states sufficient
% statistics as described in
% "A bayesian foundation for individual learning under uncertainty"
% IN:
%   - x: the previous posterior sufficient statistics. These include, using
%   the notation of [Mathys et al. 2010]:
%   mu1 = x(1);
%   mu2 = x(2);
%   sa2 = x(3);
%   mu3 = x(4);
%   sa3 = x(5);
%   Note: x(6) and x(7) are dummy states that serve in emitting a response
%   at each trial (c.f. response model) and do not play a role in the VB
%   update rules.
%       fx(6) = x(2);
%       fx(7) = x(3);
%   - P: the perceptual model parameters vector, ie. P = [ka;om;th], using
%   the notation of [Mathys et al. 2010].
%   - u: the current input to the learner.
%   - in: options set in options.inF
% OUT:
%   - fx: the updated posterior sufficient statistics (having observed u),
%   as well as book keeping of the current belief (c.f. response model).


% General remarks on notation
% - pe : prediction error
% - s : sigma
% - ka : kappa
% - om : omega
% - Xh : X 'hat'

% transform and define states and parameters
x(3) = exp(x(3)); % sa2 : sigma2 needs to be positive -> exponential mapping
x(5) = exp(x(5)); % sa3 : sigma3 needs to be positive -> exponential mapping

ka = sgm(P(1),in.kaub); % ka : kappa sigmoidal mapping + upper bound
om = P(2); % om : omega : no mapping
th = sgm(P(3),in.thub); % th : theta sigmoidal mapping + upper bound

rf = in.rf;

fx = zeros(size(x));

% 1st level
fx(1) = u(1); % [Eq. 21]

% 2nd level
s1h = sgm(x(2),1)*(1-sgm(x(2),1)); % [Eq. 26]
pe1 = fx(1) - sgm(x(2),1); % [Eq. 25] (through 21,24)
s2h = x(3) + exp(ka*x(4)+om); % [Eq. 27]
fx(3) = 1/(s2h^-1 + s1h); % [Eq. 22]
fx(2) = x(2) +fx(3)*pe1; % [Eq. 23]

% 3rd level
pi3h = 1/(x(5)+th); % [Eq. 31]
w2 = exp(ka*x(4)+om)/(x(3)+exp(ka*x(4)+om)); % [Eq. 32]
r2 = (exp(ka*x(4)+om)-x(3))/(exp(ka*x(4)+om)+x(3)); % [Eq. 33]
pe2 = (fx(3)+(fx(2)-x(2))^2)/(exp(ka*x(4)+om)+x(3)) -1; % [Eq. 34]
fx(5) = 1/(pi3h + .5*ka^2*w2*(w2+r2*pe2)); % [Eq. 29]

% Extra treatment for the variance at level 3.
% This variance can't get too small.
% When update lead to a negative variance
% this update is discarded.
% Variance is reduced (as it also is with the discarded update) through a
% rescale. This is the Rescale Factor : rf
if fx(5) <= 0
    if rf <= 0
        fx(4:5) = NaN;
    else
        fx(5) = 1/rf*x(5);
        fx(4) = x(4) + .5*ka*w2*fx(5)*pe2; % [Eq. 30]
    end
else
        fx(4) = x(4) + .5*ka*w2*fx(5)*pe2; % [Eq. 30]
end

% retransform states
x(3)  = log(x(3));
fx(3) = log(fx(3));
fx(5) = log(fx(5));

% keep track of previous reprensatation for prediction purposes
fx(6) = x(2);
fx(7) = x(3);


% TO DO:
% - include gradients dfdth, dfdx (TAKE TRANSPOSE: f is a row vector, -> option.checkGrads=1
%   for debugging)
