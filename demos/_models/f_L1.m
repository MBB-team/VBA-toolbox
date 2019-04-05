function [fx,dfdx,dfdP] = f_L1(x,P,u,in)
fx = (1+in.dt*P(1))*x + in.dt*P(2)*u;
dfdx = 1+in.dt*P;
dfdP = in.dt*[x,u];