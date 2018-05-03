function [fx] = f_FitzHughNagumo_calcium(Xt,Theta,ut,inF)

% FitzHugh-Nagumo membrane potential evolution function (+convolution)
% function [fx] = FitzHughNagumo(Xt,Theta,ut,inF)
% IN:
%   - Xt: system's states, ie:
%       Xt(1): calcium imaging prediction
%       Xt(2): membrane depolarization (mV)
%       Xt(3): proxy for ion channel opening probabilities
%   - Theta: evolution parameters (see bellow)
%   - ut: input current
%   - inF: [optional]
% OUT:
%   - fx: the evolution function evaluated at Xt

deltat = inF.delta_t;

try
    a = inF.a.^-1;
catch
    a = 1;
end
a = a.*exp(Theta(1));
K1 = (1/3)*exp(Theta(2));
K2 = 0.08*exp(Theta(3));
K3 = 0.7*exp(Theta(4));
K4 = 0.8*exp(Theta(5));
tmp = [ -a*Xt(1)+ exp(Xt(2))
        Xt(2) - K1*Xt(2).^3 - Xt(3) + ut
        K2*(Xt(2) - K3 - K4*Xt(3))];

fx = Xt + deltat.*tmp;

% A = [ -a  exp(Xt(2))
%        0  -b          ];
% dF_dX = eye(2) + deltat.*A';
% 
% dF_dTheta = [ -deltat.*Xt(1)*a  0
%               0                 -deltat.*Xt(2)*b];