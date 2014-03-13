function x = mapODT(x)
 % [useless]
a = 1;
df = -1;
cdot = exp(x);
x = (log(-a*df) + 2*log(cdot))./a;