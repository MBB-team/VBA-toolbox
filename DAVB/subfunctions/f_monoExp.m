function [fx] = f_monoExp(Xt,Theta,ut,inF)
% Hodgkin-Hoxley membrane potential evolution function (+convolution)
% function [fx] = f_calcium(Xt,Theta,ut,inF)
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
a = exp(Theta(1));
xdot = -a*Xt(1) + ut;
fx = Xt + deltat.*xdot;

