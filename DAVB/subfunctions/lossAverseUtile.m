function lau = lossAverseUtile(x,lambda)
% nu = lambda.*ones(size(x));
% nu(x>0) = 1-lambda;
% lau = 1./(1+exp(-nu.*x));

if numel(lambda)==1
    P = lambda;
    lambda(1) = 2*P;
    lambda(2) = 2*(1-P);
end
lau = zeros(size(x));
lau(x>0) = lambda(1)*log(1+x(x>0));
lau(x<0) = -lambda(2)*log(1-x(x<0));
lau(x==0) = 0;