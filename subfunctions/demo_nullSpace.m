% this demo exemplifies the projection to the null space of a GLM

clear all
close all

n = 64; % #samples
p = 8; % #regressors
X = randn(n,p);
b = abs(randn(p,1));
y0 = X*b;
y = y0 + 1e2*randn(size(y0));

c = eye(p);
[pv,stat,df,all] = GLM_contrast(X,y,c,'F',1);

bhat = all.b;
y0hat = X*bhat;
y_null = y-y0hat;
[pv,stat,df,all] = GLM_contrast(X,y_null,c,'F',1);
