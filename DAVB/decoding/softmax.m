function y = softmax(x,beta)
    y=zeros(size(x));
    for t=1:size(x,2)
        ex = exp(beta*x(:,t));
        y(:,t)=ex/sum(ex);
    end
end