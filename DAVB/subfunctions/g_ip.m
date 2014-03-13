function gx = g_ip(x,P,u,in)

gx = P(1) + real((u+P(3)).^(1-P(2)));

if isnan(gx)
    gx
end
