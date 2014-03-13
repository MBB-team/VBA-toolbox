function f = trypersistent(x,y)

persistent t xt

if isempty(t)
    t = 1;
else
    t = t+1;
end


if ~isempty(y)
    t = y;
end

xt = x + t;
f = [xt;t];