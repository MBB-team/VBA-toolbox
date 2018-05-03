function U = util(x,P)
x = x(:)';
xp = x(x>=0);
xn = x(x<0);
Uxp = P(1)*log(xp+1);
Uxn = -P(2)*log(1-xn);
U = [Uxn,Uxp];