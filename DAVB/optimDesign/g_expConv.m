function [gx] = g_expConv(x,P,u,in)

if in.up
    gx = P(1).*(1-exp(-P(2).*in.x(:)-P(3)));
else
    gx = P(1) + exp(-P(2).*in.x(:)-P(3));
end
