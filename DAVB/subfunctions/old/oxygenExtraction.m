function [y,dydx] = oxygenExtraction(x,E0)

try
    E0;
catch
    E0 = 0.34;
end
% x = exp(x);
y = (1-(1-E0).^(1./x))./E0;
dydx = log(1-E0).*(1-E0).^(1./x)./(E0.*x);

alpha = 2;
y = x.^(1./alpha);                        % blood outflow
dydx = (1./alpha).*y;


[E0,dsdp] = sigm(x-0.6633,struct('beta',2));
E0 = E0(:);
dsdp = dsdp(:);

ff = (1-(1-E0).^(1./2))./E0;
y = (2.*ff - 3)./(3);

% dydx = 2.*(0.5*E0.*(1-E0).^(-1+1./2)-(1-(1-E0).^(1./2))).*dsdp./(3.*E0.^2);
dydx = 2.*(0.5*(1-E0).^(-1+1./2)-ff).*dsdp./(3.*E0);


% y = (1-E0).^(1./2);
% dydx = -0.5.*(1-E0).^(-1+0.5).*dsdp;

