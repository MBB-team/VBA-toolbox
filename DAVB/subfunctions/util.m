function U = util(x)
x = x(:)';
xp = x(x>=0);
xn = x(x<0);
Uxp = log(xp+1);
Uxn = -log(1-xn);
U = [Uxn,Uxp];