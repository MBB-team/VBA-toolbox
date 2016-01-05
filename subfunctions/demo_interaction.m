% this demo studies the non-trivial properties of statistical interactions
% Let us consider the following scenario:
% - we construct an experimental design with two conditions
% - we measure two variables in both conditions: our dependant measure (y)
% and another (potentially modrating) variables X1.
% Now we want to know whether different analyses may help us in identifying
% whether X1 modulates the efefct of the experimental manipulation onto y.

clc
clear all
close all

n = 32;
X1 = [ones(n,1);-ones(n,1)];
X2 = randn(2*n,1);
Xi = X1.*X2;
X = [X1,X2,Xi,ones(2*n,1)];
b = [1;1;10;1];
e = 0.1*randn(2*n,1);
y = X*b + e;

[p] = GLM_contrast(X,y,[0;0;1;0],'F',1);

% now let us look at X2 in terms of a contrast:
yc = y(X1==1) - y(X1==-1);
x1c = X2(X1==1) - X2(X1==-1);
[p] = GLM_contrast([x1c,ones(n,1)],yc,[1;0],'F',1);

