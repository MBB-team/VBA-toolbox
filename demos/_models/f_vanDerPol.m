function [fx,dF_dX,dF_dTheta] = f_vanDerPol(Xt,Theta,ut,inF)
% evolution function for the Van Der Pol (limit cycle) oscillator.
% function [fx,dF_dX,dF_dTheta,d2F_dXdTheta] = f_vanDerPol(Xt,Theta,ut,inF)
% IN:
%   - Xt: current hidden state
%   - Theta: evolution parameter
%   - ut: current input to the system
%   - in: user-defined input structure (containing the time decimation)
% OUT:
%   - fx: the predicted hidden state
%   - dF_dX: the jacobian of the system
%   - dF_dTheta: the derivative of the evolution function w.r.t the
%   evolution parameters

deltat = inF.deltat;

x       = Xt;
mu       = Theta(1);

fx      = zeros(2,1);
fx(1)   = x(2);
fx(2)   = mu.*(1-x(1).^2).*x(2) - x(1);
fx      = deltat.*fx + x;

J       = zeros(2,2);
J(1,:)  = [0,1];
J(2,:)  = [-2.*mu.*x(1).*x(2)-1,mu.*(1-x(1).^2)];
dF_dX   = deltat.*J';

dF_dTheta           = zeros(1,2);
dF_dTheta(1,:)      = [0,(1-x(1).^2).*x(2)];
dF_dTheta           = deltat.*dF_dTheta;

                       