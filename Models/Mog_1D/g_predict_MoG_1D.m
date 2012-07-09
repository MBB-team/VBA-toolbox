function [gx] = g_predict_MoG_1D(x,P,u,in)



% alpha = sigm(P(1));
% gam = exp(P(2));
% mu0 = P(3);
% s0 = exp(P(4));

pi1 = sigm(P(1));
mu = [P(2),P(4)];
s = exp([P(3),P(5)]);
y = u(1);

c = gp(y,mu(1),s(1))/gp(y,mu(2),s(2))*pi1/(1-pi1);
p1 = c/(1+c); % prob class 1


gx = gp(y,mu(1),s(1))*pi1 + gp(y,mu(2),s(2))*(1-pi1)

   

function gp = gp(x,mu,s)
gp = 1/sqrt(2*pi)/s*exp(-1/2*(x-mu)^2/s^2);
