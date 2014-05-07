function [gx] = g_dummy(x,P,u,in)

gx = in.X*P(1);
try
    gx = gx + in.X*P(2);
end
