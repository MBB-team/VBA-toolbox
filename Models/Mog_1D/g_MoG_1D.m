function [gx] = g_MoG_1D(x,P,u,in)

pi1 = sigm(P(1));

mu = [P(2),P(4)];
s = exp([P(3),P(5)]);
i = (rand<pi1) + 1;
gx = mu(i) + randn*s(i);
    
