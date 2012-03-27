function [r1,r2,r3,e1,e2] = designEfficiency_GLM(X,c,alpha,sigma)
% checks consistency of design efficiency expressions
% [r1,r2,r3] = designEfficiency_GLM(X,c,alpha,sigma)
% IN:
%   - X: GLM design matrix
%   - c: contrast vector (assumed here to be all-weros but one entry)
%   - alpha: prior precision of model parameters
%   - sigma: residuals precision
% OUT:
%   - r1,r2,r3: three expressions for the design Bayesian risk.
%   Note: these expressions should give the exact same result.
%   - e1,e2: two expressions for the clqssicql design efficiency.
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

n = size(X,1);
X1 = X;
X2 = X(:,c==0);
Xi = X(:,c==1);

Q1 = eye(n)/sigma + X1*X1'/alpha;
Q2 = eye(n)/sigma + X2*X2'/alpha;


a = Xi'*inv(Q2)*Xi./alpha;
r1 = -0.25*log(1+(a.^2./(4*(1+a))));

r2 = -0.5*log(det(0.5*(Q1+Q2))) ...
    +0.25*log(det(Q1)) + +0.25*log(det(Q2));

r3 = -0.5*log(det(Q2+0.5* Xi*Xi'/alpha)...
    ./sqrt(det((Q2+ Xi*Xi'/alpha)*Q2)));

e1 = 1./(c'*inv(X'*X)*c*sigma);

e2 = (Xi'*(eye(n)-X2*inv(X2'*X2)*X2')*Xi)./sigma;


