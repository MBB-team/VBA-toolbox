function gx = softmaxMat(x)
% sotmax mapping, i.e.: g(x) = exp(x)./sum(exp(x))
% function gx = softmax(x)
x = x-repmat(max(x,[],1),size(x,1),1);
gx = exp(x)./repmat(sum(exp(x),1),size(x,1),1);