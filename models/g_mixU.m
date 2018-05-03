function [gx] = g_mixU(x,P,u,in)
P = P.*in.weights;
gx = P(1)*u(1) + P(2)*(u(1)+u(2)) + u(3)*(P(3)+P(4)) + u(4)*P(5) + P(6);