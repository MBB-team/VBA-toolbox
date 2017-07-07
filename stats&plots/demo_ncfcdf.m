
clear all
close all
clc
for i=1:Nmc
    x = randn.^2;
    nu1 = floor(1+64*rand);
    nu2 = nu1+floor(1+64*rand)+2;
    p(i,1) = ncfcdf(x,nu1,nu2,delta);
    p(i,2) = VBA_ncfcdf(x,nu1,nu2,delta);
end

figure,plot(p(:,1),p(:,2),'.')
