function  [ fx ] = f_VBfree( x,P,u,in )
% IN:
% - x_t : two posterior moments of the mean and variance of u^(o)
% - P : volatility of the moments (2)
% - u_t : previous feedback
% - in : []

theta = exp(P);

o = u(2);
ia = 4*(1-u(1));

mu1 = x(1+ia);
s1 = x(2+ia);
mu2 = x(3+ia);
s2 = x(4+ia);

Eexpx2 = exp(-mu2);%+s2/2);
s1 = 1./(Eexpx2+1./(s1+theta(1)));
mu1 = mu1 + s1.*Eexpx2*(o-mu1);

Esquarederr = (o-mu1).^2 + s1;
s2 = 1./(exp(-mu2)*Esquarederr/2+1./(s2+theta(2)));
mu2 = mu2 - s2.*(1-exp(-mu2)*Esquarederr)/2;

fx = x;
fx(ia+[1:4]') = [mu1;s1;mu2;s2];