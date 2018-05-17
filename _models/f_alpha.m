function [fx] = f_alpha(x,P,u,in)
H = P(1);
T = exp(P(2));

dx = [x(2); ...
      (H / T) * u - (2 / T) * x(2) - ( 1 / T ^ 2) * x(1)];
  
fx = x + in.dt * dx;