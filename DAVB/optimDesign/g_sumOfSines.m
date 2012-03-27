function [gx,dgdx,dgdP] = g_sumOfSines(x,P,u,in)

ts = in.t0:in.dt:in.tf;
X = [sin(ts(:)),sin(ts(:)+u(1))];
gx = X*P;
dgdx =[];
dgdP = X';