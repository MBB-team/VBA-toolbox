function [fx] = f_HH_calcium(Xt,Theta,ut,inF)
% Hodgkin-Hoxley membrane potential evolution function (+convolution)
% function [fx] = f_HH_calcium(Xt,Theta,ut,inF)
% IN:
%   - Xt: system's states, ie:
%       Xt(1): calcium imaging prediction
%       Xt(2): membrane depolarization (mV)
%       Xt(3:5): type I-II-III ion channel opening probabilities, in
%       inverse-sigmoid space (Gaussian transformation)
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
C = exp(Theta(2));
gK = 36*exp(Theta(3));
gNa = 120*exp(Theta(4));
gL = 0.3*exp(Theta(5));

EK = -12;
ENa = 115;
EL = 10.6;

V = Xt(2);
m = sigm(Xt(3));
n = sigm(Xt(4));
h = sigm(Xt(5));

an = (0.1 - 0.01*V)./(exp(1-0.1*V)-1);
am = (2.5 - 0.1*V)./(exp(2.5-0.1*V)-1);
ah = 0.07*exp(-V./20);
bn = 0.125*exp(-V/80);
bm = 4*exp(-V/18);
bh = 1./(exp(3-0.1*V)+1);

tmp = [ -a*Xt(1)+ V*1e-0
        (-gNa*m.^3*h*(V-ENa) - gK*n.^4*(V-EK) - gL*(V-EL) + 1e0*ut)/C
        (am*(1-m) - bm*m)/(m-m.^2+1e-2)
        (an*(1-n) - bn*n)/(n-n.^2+1e-2)
        (ah*(1-h) - bh*h)/(h-h.^2+1e-2)   ];

fx = Xt + deltat.*tmp;

