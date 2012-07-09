function [gx] = g_classify_MoG_1D(x,P,u,in)


pi1 = sigm(P(1));
mu = [P(2),P(4)];
s = exp([P(3),P(5)]);

y = u(1);

c = gp(y,mu(1),s(1))/gp(y,mu(2),s(2))*pi1/(1-pi1);
%gx = c>1;

gx = c/(1+c);
    

function gp = gp(x,mu,s)
gp = 1/sqrt(2*pi)/s*exp(-1/2*(x-mu)^2/s^2);
