function [fx,dfdx] = f_LV2D(x,P,u,in)


a = P(1);
b = P(2);
c = P(3);
d = P(4);

dxdt = [a*x(1) - b*x(1)*x(2);-c*x(2) + d*x(1)*x(2)];

fx = x + in.deltat*dxdt;
dfdx = [    a-b.*x(2)   -b*x(1)
            d*x(2)      d*x(1)-c    ];
dfdx = eye(length(x)) + in.deltat*dfdx;
