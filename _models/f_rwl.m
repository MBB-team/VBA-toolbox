function fx = f_rwl(x,P,u,in)
fx = x;
for i=1:size(u,1)
    fx = fx + P(i)*(u(i)-x);
end