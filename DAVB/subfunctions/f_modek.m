function [fx,dfdx] = f_modek(x,P,u,in)
% evolution function for dynamical modes of neural field

if in.k == 0
    z_10 = x(1);
else
    z_10 = in.z_10;
end

v = in.v(in.i1,in.i2);
c = in.c(in.i1,in.i2);
a = in.a(in.i1,in.i2);

cz = 0.5*dsdv(z_10,in).*v.*a./c;
nk = length(in.k);
J = zeros(nk*3,nk*3);
for k=1:nk
    
    K = -(in.k(k)*pi/in.l)^2;
    fk = 0.5*v*(K*c^2-1)/c;
    
    Jk = [  0           0               1
            cz          fk              0
            -in.Ke^2    in.me*in.Ke     -2*in.Ke  ];
    
    J((k-1)*3+1:k*3,(k-1)*3+1:k*3) = Jk;
    
end

% hf=figure,imagesc(J)
% pause
% close(hf)

xdot = J*x + in.C*u;

fx = x + in.dt*xdot;
dfdx = eye(3*nk) + in.dt*J;