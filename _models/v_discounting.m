function [v] = v_discounting(x,P,u,in)
% evaulates 2D utility function at u
% function [v] = v_discounting(x,P,u,in)
% This function computes 2D utility functions, as in, e.g., temporal
% discounting models. Four different functional forms can be used:
% hyperbolic, exponential, linear and basis functions sets.
% IN:
%   - x: [useless]
%   - P: discounting parameter
%   - u: 2x1 vector of utility features (e.g., in the context of temporal
%   discounting, u(1)=time, u(2)=reward).
%   - in: structure containing the following fields:
%       .model: can be 'hyperbolic', 'exponential', 'linear' or 'basis'
%       .ind.logk: index of the discount parameter (for the three first
%       model types)
%       .gx and .gy [only for .model='basis']: grids for 1st and 2nd
%       dimensions of the utility functions
%       .bf [only for .model='basis']: 3D array contining the basis
%       functions
% OUT:
%   - v: utility function, evaluated at u.

% define dimensions
t = u(1);
R = u(2);

% evaluates utility function at u
switch in.model
    case 'hyperbolic'
        k = exp(P(in.ind.logk));
        v = R./(1+k.*t);
    case 'exponential'
        k = exp(P(in.ind.logk));
        v = R.*exp(-k.*t);
    case 'linear'
        k = P(in.ind.logk);
        v = R - k*t;
    case 'basis'
        [tmp,i1] = min((in.gx-t).^2);
        [tmp,i2] = min((in.gy-R).^2);
        v = VBA_vec(in.bf(i1,i2,:))'*P;
end

