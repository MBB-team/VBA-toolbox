function y=sig(x)
y = 1./(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;