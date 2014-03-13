function [I] = dummyf(x)
% dummy bimodal density
I = -(1/8).*(x(1)-(1/16)*x(2).^2).^2 - (1/128)*(x(2).^2-2).^2;
I = exp(I);