function [y,dy] = safePositif(x)
% safePositif: max(0,x) but toolbox friendly (with gradient)
% y      = safePositif(x)
% [y,dy] = safePositif(x)

k = 150;

%%
y = log(1+exp(k*x))/k;
y(isinf(y))=x(isinf(y));

%%
if nargout>1
    dy=exp(k*x) ./ (exp(k*x) + 1); 
    dy(isnan(dy))=1;
end

end
