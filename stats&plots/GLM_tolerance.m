function tol = GLM_tolerance(X0)
% computes regression tolerance of a given design matrix
% function tol = GLM_tolerance(X)
% In the context of multiple regression, tolerance measures how well each
% regression weight can be estimated, given that the design matrix may
% suffer from multicollinearity. Note: multicollinearity does not reduce
% the reliability of the model as a whole, it only affects calculations
% regarding individual predictors. That is, a multiple regression model
% with correlated predictors can indicate how well the entire bundle of
% predictors predicts the outcome variable, but it may not give valid
% results about any individual predictor.
% IN:
%   - X0: nxp design matrix, where n is the number of data samples and p
%   is the number of regressors.
% OUT:
%   - tol: px1 tolerance vector. By convention, a tolerance of less than
%   0.10 indicates a multicollinearity problem.

[n,p] = size(X0);

% z-score
X0 = bsxfun(@minus,X0, mean(X0));
sigma0 = std(X0);
sigma0(sigma0==0) = 1;
X0 = bsxfun(@rdivide,X0, sigma0);

tol = zeros(p,1);
for i=1:p
    X = X0(:,setdiff(1:p,i));
    y = X0(:,i);
    C = X'*X;
    iC = pinv(C);
    b = iC*X'*y;
    yhat = X*b;
    SS_tot = sum((y-mean(y)).^2);
    SS_err = sum((y-yhat).^2);
    tol(i) = SS_err/SS_tot;
end