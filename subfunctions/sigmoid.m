function   y = sigmoid(x)
y=1./(1+exp(-x));
if y < 1e-4
    y = 1e-4*ones(size(y));
elseif y > 1-1e-4
    y =( 1-1e-4)*ones(size(y));
end