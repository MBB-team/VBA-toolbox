function [gx] = g_GammaDensity(x,P,u,in)
gamma = exp(P(1));
tau = exp(P(2));
K = P(3);
gx = K + gamma.*tau.*in.grid.*exp(-tau*in.grid);
