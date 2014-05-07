
close all
clear all

% simulate dummy independant and dependant variables in condition 1
n1 = 32;
x1 = randn(n1,1);
y1 = x1 + randn(n1,1);


% simulate dummy independant and dependant variables in condition 2
n2 = 32;
x2 = randn(n2,1);
y2 = 0.2*x2 + randn(n2,1);

% form GLM
Y = [normalize(y1(:));normalize(y2(:))];
X = zeros(n1+n2,2);
X(1:n1,1) = normalize(x1(:));
X(n1+1:n1+n2,2) = normalize(x2(:));

% test signed difference
c = [1;-1];
type = 't';
verbose = 1;
Xnames = {'cond 1','cond 2'};
[pv,stat,df,all] = GLM_contrast(X,Y,c,type,verbose,Xnames)

% test unsigned difference
c = [1;-1];
type = 'F';
verbose = 1;
Xnames = {'cond 1','cond 2'};
[pv,stat,df,all] = GLM_contrast(X,Y,c,type,verbose,Xnames)