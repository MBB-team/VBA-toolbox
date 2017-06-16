function sse = SSE(x,y)
% sum-of-squared distance between x and y
sse = sum((vec(x)-vec(y)).^2);