function [fx,dF_dX,dF_dTheta,d2F_dXdTheta] = f_doubleWell(Xt,Theta,ut,inF)
% damped double Well Lyapunov evolution function
% function [f,J] = f_doubleWell(t,x,theta)
%
% This function computes the evolution function that comes from a .
% IN:
%   - t: time index (not used here)
%   - x: the current state of the system
%   - theta: a 2x1 vector parameter (containing the position of the two
%   wells.
% OUT:
%   - f: the current value of the evolution function
%   - J: the jacobbian of the system
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

deltat = inF.deltat;

a   = Theta(1);
b   = Theta(2);
k   = Theta(3);
x1  = Xt(1);
x2  = Xt(2);


f   = [ x2  ;  -2*(x1-a)*(x1-b)^2-2*(x1-a)^2*(x1-b)-k.*x2 ];
J   = [ 0                                               1
        -2*(x1-b)^2-8*(x1-a)*(x1-b)-2*(x1-a)^2        -k ];


fx                  = deltat.*f + Xt;
dF_dX               = deltat.*J' + eye(2);
dF_dTheta           = deltat.*[ 0 2*(x1-b)^2+4*(x1-a)*(x1-b)
                                0 4*(x1-a)*(x1-b)+2*(x1-a)^2
                                0 -x2                        ];
d2F_dXdTheta(:,:,1) = zeros(2,3);
d2F_dXdTheta(:,:,2) = deltat.*[ 12*x1-8*b-4*a   12*x1-4*b-8*a   0 
                                0               0               -1 ];



