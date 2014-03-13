function x = mapDT(x,theta)
% [useless]
alpha = theta(1);
dF = theta(2);
a = theta(3);

x = (log(alpha.*dF) - 2*log(a*x))./alpha;
