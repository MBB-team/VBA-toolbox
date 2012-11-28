function [gx,dgdx,dgdP] = g_x2(x,P,u,in)

try x=in.x(:); catch; x=1; end

gx = 1 + 0.1*x.*P - x.*P.^2;% + x.*P.^3;
dgdx = [];
dgdP = 0.1*x' - 2*x'.*P;% + x'.*3.*P.^2;