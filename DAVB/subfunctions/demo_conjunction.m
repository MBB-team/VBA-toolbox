% This demo exemplifies a conjunction analysis, i.e. a joint test that two
% parameters are both either positive or negative.

clear all
close all
n = 1e2;


X1 = randn(n,1); % reward
X2 = randn(n,1); % effort
X = [X1,X2,ones(n,1)];


b1 = [-1;1;randn];
e = randn(n,1);
y1 = X*b1 + e;

b2 = [1;1;randn];
e = randn(n,1);
y2 = X*b2 + e;

y = [y1,y2];
y = repmat(y,1,4);

c = eye(3);
[pv,stat,df,all] = GLM_contrast(X,y,c,'F',1,{'reward','effort','mean'});

C = {[1;0;0],[0;1;0],[-1;0;0],[0;-1;0]};

for i=1:4
    [p1(i),stat,df,all] = GLM_contrast(X,y1,C{i},'t',0);
    [p2(i),stat,df,all] = GLM_contrast(X,y2,C{i},'t',0);
end

% form conjunction, ie test whether reward and effort effects go in the
% same direction (for both datasets)
pc1 = 1- ( (1-p1(1))*(1-p1(2))+(1-p1(3))*(1-p1(4)) )
pc2 = 1- ( (1-p2(1))*(1-p2(2))+(1-p2(3))*(1-p2(4)) )


