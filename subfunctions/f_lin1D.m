function [fx,dF_dX,dF_dTheta] = f_lin1D(Xt,Theta,ut,inF)
% dummy 2D linear evolution function

deltat = inF.delta_t;

try
    a = inF.a.^-1;
catch
    a = 1;
end
try
    [uu,dudtheta] = feval(inF.u_fname,Theta(2:end),ut(2:end),inF);
catch
    uu = 0;
    if size(Theta,1) > 1
        dudtheta = zeros(size(Theta,1)-1,1);
    else
        dudtheta = [];
    end
end

A = -a.*exp(Theta(1));
fx = Xt + deltat.*(A*Xt + ut(1) + uu);
dF_dX = 1 + deltat.*A;
dF_dTheta = deltat.*([Xt*A;zeros(size(Theta,1)-1,1)] + [0;dudtheta]);
