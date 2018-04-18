function [E, V] = VBA_dirichlet_moments (a)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [E, V] = VBA_dirichlet_moments (a)
% Derives the first- and second-order moments of a Dirichlet density
%
% IN:
%   - a: parameters of the Dirichlet distribution (pseudo-counts)
% OUT:
%   - E: first-order moment (expectation)
%   - V: second-order moment (variance)
%
% /////////////////////////////////////////////////////////////////////////

% Let X = (X_1, ..., X_k) ~ Dir(a)

% shorthand
a0 = sum(a);

% First-order moment E[X]
E = a./a0;

% second order moment Var[X]
V = (a .* (a0-a)) ./ (a0.^2 * (a0+1)); 