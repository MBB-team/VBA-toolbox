function [gx,dgdx,dgdp] = g_AVL(x,Phi,u,in)
% OTO: reaction time observation function
% function [gx,dgdx,dgdp] = g_AVL(x,Phi,u,in)
% This function computes the observation function (reaction time) for the
% audio-visual associative learning task.
% IN:
%   - x: the 5x1 vector of states. The last one is the only relevant state
%   for the response model: it is the 'VB prediction error', i.e. the
%   difference between the posterior and the prior 1st-order moments of the
%   cue identity (the smaller it is, the quicker the reaction time).
%   - Phi: 2x1 vector of observation parameters. Phi(1) is an intercept
%   constant and Phi(2) is a scaling factor.
%   - u: of importance here is the choice of the subect at the curren
%   trial, since a categorization error leads to a quick reaction time.
%   - in: contains the index of the subject's choice in the u vector.
% OUT:
%   - gx: expected reaction time
%   - dgdx: derivative wrt states

y = 2*exp(Phi(1)-Phi(2)).*(x(1)-x(end)).*(2*u(in.uc)-1);
le = y;
dledy = 1;

if le > 0
    
    gx = 0.5.*exp(-Phi(1)).*log(le);
    
    dydx = [    2*exp(Phi(1)-Phi(2)).*(2*u(in.uc)-1)
        zeros(size(x,1)-2,1)
        -2*exp(Phi(1)-Phi(2)).*(2*u(in.uc)-1) ];
    dgdx = 0.5.*exp(-Phi(1)).*dledy.*dydx./le;
    
    dydp = [ 2*exp(Phi(1)-Phi(2)).*(x(1)-x(end)).*(2*u(in.uc)-1)
        -2*exp(Phi(1)-Phi(2)).*(x(1)-x(end)).*(2*u(in.uc)-1) ];
    dgdp = [ 0.5.*exp(-Phi(1)).*(-log(le) + dledy.*dydp(1)./le)
        0.5.*exp(-Phi(1)).*dledy.*dydp(2)./le              ];

    if size(Phi,1) > 2
        gx = gx + Phi(3);
        dgdp = [ dgdp ; 1 ];
    end
    
else
    
    gx = 0;
    dgdx = zeros(size(x));
    dgdp = zeros(size(Phi));
    
end


