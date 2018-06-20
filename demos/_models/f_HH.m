function [fx] = f_HH(Xt, Theta, ut, inF)
% Hodgkin-Hoxley membrane potential evolution function
% function [fx] = f_HH(Xt,Theta,ut,inF)
% IN:
%   - Xt: system's states, ie:
%       Xt(1): membrane depolarization (mV)
%       Xt(2:4): type I-II-III ion channel opening probabilities, in
%       inverse-sigmoid space (Gaussian transformation)
%   - Theta: evolution parameters (see bellow)
%   - ut: input current
%   - inF: mandatory for field .delta_t for time discretization
% OUT:
%   - fx: the evolution function evaluated at Xt


deltat = inF.delta_t;
C = exp(Theta(1)); % neuron capacitance
gK = 36*exp(Theta(2));
gNa = 120*exp(Theta(3));
gL = 0.3*exp(Theta(4));

EK = -12;
ENa = 115;
EL = 10.6;

V = Xt(1);
m = VBA_sigmoid(Xt(2));
n = VBA_sigmoid(Xt(3));
h = VBA_sigmoid(Xt(4));

an = (0.1 - 0.01*V)./(exp(1-0.1*V)-1);
am = (2.5 - 0.1*V)./(exp(2.5-0.1*V)-1);
ah = 0.07*exp(-V./20);
bn = 0.125*exp(-V/80);
bm = 4*exp(-V/18);
bh = 1./(exp(3-0.1*V)+1);

r = 1e-2; % for numerical stability

xdot = [ (-gNa*m.^3*h*(V-ENa) - gK*n.^4*(V-EK) - gL*(V-EL) + ut)/C
        (am*(1-m) - bm*m)/(m-m.^2+r)
        (an*(1-n) - bn*n)/(n-n.^2+r)
        (ah*(1-h) - bh*h)/(h-h.^2+r)   ];

fx = Xt + deltat.*xdot;

end
