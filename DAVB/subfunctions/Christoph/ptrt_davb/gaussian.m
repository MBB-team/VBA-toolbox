% -------------------------------------------------------------------------
function p = gaussian(x,mu,sigma)
    p = 1/sqrt(2.*pi.*sigma).*exp(-(x-mu).^2./(2.*sigma));
return;
