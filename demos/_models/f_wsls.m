function  [fx] = f_wsls(x,P,u,in)
% pseudo-evolution function for a "win/stay - lose/switch" strategy
% function  [fx,dfdx,dfdP] = f_wsls(x,P,u,in)
% The "win/stay - lose/switch" strategy is encoded in terms of the
% evolution of pseudo q-values, which swap sign depending upon the
% feedback. Note: the feedback should unambiguously code for "win" or
% "lose", i.e. in terms of a positive (resp. negative) reward.
% IN:
%   - x : pseudo q-values (1: stay, -1:switch)
%   - P : [useless]
%   - u : u(1)=previous action, u(2)=feedback
%   - in : [useless]
% OUT:
%   - fx: evolved pseudo q-values (2x1)
%   - dfdx/dfdP: [useless]

r = u(2);
if u(1)==1
    fx = sign(r)*[1;-1];
else
    fx = sign(r)*[-1;1];
end
