function gx = softmax(x,beta)
% sotmax mapping, i.e.: g(x) = exp(x)./sum(exp(x))
% function gx = softmax(x,beta)
try,beta;catch,beta=1;end
x = x-max(x);
gx = exp(beta*x)./(sum(exp(beta*x)));