function [y,labels] = generateGMM(n,lambda,mu,gamma,verbose)
%
% function y = generateGMM(n,lambda,mu)
% This function generates n observations following a MoG distribution
% given by the class frequencies lambda, their mean mu and their variance gamma. 

try, verbose = ~~verbose; catch, verbose = 1; end

p = size(mu,1);
k = 1;
sumLambdak = lambda(1);

if verbose && (p == 2 || p == 3)
    col = repmat('bgrcmyk',1,ceil(1+length(gamma)./7));
    pl = figure;
    hold on
end

labels = zeros(n,length(gamma));

for i = 1:n
    
    ratio = i./n;
    
    if ratio>=sumLambdak && ~isequal(i,n)
        k = k+1;
        sumLambdak = sumLambdak+lambda(k);
    end
    
    labels(i,k) = 1;
    y(:,i) = mu(:,k) + sqrt(gamma(k)).*randn(p,1);
    
    if verbose && p == 2
        figure(pl),plot(y(1,i),y(2,i),[col(k),'.']);
    elseif verbose && p == 3
        figure(pl),plot3(y(1,i),y(2,i),y(3,i),[col(k),'.']);
    end
%     drawnow;
end


if verbose &&(p==2 || p==3)
    axis square
end