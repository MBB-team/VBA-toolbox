function [ gx ] = g_lm( x,P,u,in )
%
% u :inputs
% P :modulators


gx = P(1).*u(1)+P(2).*u(2)+P(3).*u(3);

end

