function [m,v] = VBA_EVprodX(mu,Sigma)
% calculates the first two moments of the product of two random variables

% function [m,v] = EVprodX(mu,Sigma)
% IN:
%   - mu: (2x1 vector) mean of the two variables
%   - Sigma: (2x2 matrix) covariance matrix of the two variables
% OUT:
%   - m: expectation of the product of the two variables
%   - v: variance of the product of the two variables


A = [ 0 1
      1 0 ];
m = mu(1)*mu(2) + Sigma(1,2);
v = 2*trace(A*Sigma*A*Sigma) + 4*mu'*A*Sigma*A*mu;