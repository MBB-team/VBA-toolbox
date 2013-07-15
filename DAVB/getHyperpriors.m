function [a, b]=getHyperpriors(variance,p_min,p_max)
%% compute hyperpriors parameters given you want to eplain between p_min 
%  and p_max fraction of the total variance.
%% ex.
% [a, b]=getHyperpriors(200,.1,.9);
% figure;
% r=random('gam',a,1/b,1,5000);
% hist(r,50);
% hold on
% plot([prec_min, prec_max],[100 100],'r');

%% check parameters
if nargin==1
    p_min=.1;
    p_max=.9;
elseif nargin==3
    if any([p_min, p_max]<0) || any([p_min, p_max]>1)
        error('*** getHyperpriors: p_min and p_max should be between 0 and 1.');
    end
    if p_min > p_max
        error('***  getHyperpriors: p_min should be inferior to p_max.');
    end
else
    error('***  getHyperpriors: wrong number of argument');
end

p_max = min(p_max, 1-eps);

% residual variance
var_min = variance*(1-p_max);
var_max = variance*(1-p_min);

% expressed as precision
prec_min = 1/var_max;
prec_max = 1/var_min;

% approximate 98% confidence interval
a = 6*((prec_min+prec_max)^2)/((prec_max-prec_min)^2);
b = 12*(prec_min+prec_max)/((prec_max-prec_min)^2);

