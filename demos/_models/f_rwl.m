function fx = f_rwl(x,P,u,in)
fx = x;
for i = 1 : size (u, 1)
    if ~ VBA_isWeird (u(i))
        fx = fx + P(i) * (u(i) - x);
    end
end