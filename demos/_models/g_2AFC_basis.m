function [gx,dgdx,dgdp] = g_2AFC_basis(x,P,u,in)
% 2AFC observation function: projection onto basis function set
% function [gx,dgdx,dgdp] = g_2AFC_basis(x,P,u,in)
% This function evaluates the probability of choosing the first alternative
% options in a 2 alternative forced choice paradigm, where each option is
% decribed in terms of two orthogonal dimensions (e.g. reward & cost).
% Critically, the code uses a semi-parametric utility model, i.e. the
% 2D-utility plane is described using a basis function set provided in the
% 'in' input structure.
% IN:
%   - x: [useless]
%   - P: projection parameters onto the basis function set
%   - in: structure containing the following fields:
%       .gx: grid of the 1st-dim along which the utility basis functions
%       are evaluated
%       .gx: grid of the 2nd-dim
%       .ind.x/ind.y: indices of the 1st (resp., 2nd) dimension of the two
%       alternative options.
%       .bf: 3D-cell array containing the utility basis functions (one per
%       entry along the third dimension).
% OUT:
%   - gx: the probability of picking the firt alternative
%   - dgdx: [useless]
%   - dgdp: gradient of the probability wrt projection parameters

x = u(in.ind.x);
y = u(in.ind.y);

V = zeros(length(P),2);
for i=1:2
    [tmp,ix(i)] = min((in.gx-x(i)).^2);
    [tmp,iy(i)] = min((in.gy-y(i)).^2);
    V(:,i) = VBA_vec(in.bf(ix(i),iy(i),:));
end
dV = V(:,1)-V(:,2);
dv = dV'*P;
gx = VBA_sigmoid(dv);
dsdx = gx.*(1-gx);
dgdx = [];
dgdp = dsdx*dV;

