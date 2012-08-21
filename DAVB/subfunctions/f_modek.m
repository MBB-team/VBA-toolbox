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

cz = dsdv(z_10,in).*v*c.*a/2;
K = -in.k*pi/in.l;
fk = -v*c*(1-K./c^2);

J = [   0           0           1
        cz          fk          0
        -in.Ke^2    in.me*in.Ke -2*in.Ke  ];

xdot = J*x + in.C*u;

fx = x + in.dt*xdot;
dfdx = eye(3) + in.dt*J;