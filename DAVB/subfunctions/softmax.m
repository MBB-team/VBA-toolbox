function gx = softmax(x)
x = x-max(x);
gx = exp(x)./(sum(exp(x)));