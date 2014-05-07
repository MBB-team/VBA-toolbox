function [fx] = f_FitzHughNagumo(Xt,Theta,ut,inF)

% FitzHugh-Nagumo membrane potential evolution function
% function [fx] = FitzHughNagumo(Xt,Theta,ut,inF)
% IN:
%   - Xt: system's states, ie:
%       Xt(1): membrane depolarization (mV)
%       Xt(2): proxy for ion channel opening probabilities
%   - Theta: evolution parameters (see below)
%   - ut: input current
%   - inF: [optional]
% OUT:
%   - fx: the evolution function evaluated at Xt


K1 = (1/3)*exp(Theta(1));
K2 = 0.08*exp(Theta(2));
K3 = 0.7*exp(Theta(3));
K4 = 0.8*exp(Theta(4));
xdot = [    Xt(1) - K1*Xt(1).^3 - Xt(2) + ut
            K2*(Xt(1) - K3 - K4*Xt(2))      ];

fx = Xt + inF.dt.*xdot;
