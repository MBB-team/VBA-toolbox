function [fx] = f_HRF(Xt,P,ut,in)
% Balloon (HRF) model evolution function in log-space (wo options)
% function [fx] = f_HRF(Xt,P,ut,in)
% This function evaluates the evolution function derived from the balloon
% model for the hemodynamic response function. Note that the hemodynamic
% states are in log-space, for positivity constraints.


% Get parameters
[E0,V0,tau0,kaf,kas,epsilon,alpha] = BOLD_parameters;
if ~isempty(P)
    iE0 = VBA_sigmoid(E0,'inverse',true);
    E0 = VBA_sigmoid(P(1)+iE0);
    tau0 = tau0.*exp(P(2));
    kaf = kaf.*exp(P(3));
    kas = kas.*exp(P(4));
    epsilon = epsilon.*exp(P(5));
    alpha = alpha.*exp(P(6));
end

% Exponentiate hemodynamic states
x1 = Xt(1); % vasodilatory signal s(t)
x2 = exp(Xt(2)); % blood inflow f(t)
x3 = exp(Xt(3)); % blood volume v(t)
x4 = exp(Xt(4)); % dHb content q(t)

% Intermediate variables
fv = x3.^(1./alpha); % blood outflow
ff = (1-(1-E0).^(1./x2))./E0; % oxygen extraction

% Evaluate flow field
f = [   epsilon.*ut - kas.*x1 - kaf.*(x2 - 1)
        x1./x2
        (x2 - fv)./(tau0.*x3)
        (x2.*ff./x4 - fv./x3)./tau0             ];

try
    deltat = in.deltat;
    fx = Xt + deltat.*f;
catch
    fx = f;
end



