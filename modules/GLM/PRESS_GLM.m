function [R2_cv,mP,VP,R2] = PRESS_GLM(X,y,verbose)
% evaluates the PRESS-R2 cross-valuidation metric for a GLM
% function [R2_cv,mP,VP,R2] = PRESS_GLM(X,y,verbose)
% This function evaluates a prediction accuracy metric (PRESS-R2), in the
% context of a leave-one-out cross-validation scheme for a GLM. Note: it is
% based upon the predicted residual error sum of squares (PRESS) statistic.
% IN:
%   - X: nxp design matrix (independent variables)
%   - y: nxk dependent variables
%   - verbose: verbose flag
% OUT:
%   - R2_cv: PRESS-R2, i.e. percentage of predicted variance, where
%   predictions are derived for each fold of the leave-one-out
%   cross-validation scheme.
%   - mP: pxk OLS-estimates of the GLM parameters
%   - VP: pxk varianceso of the parameters estimates
%   - R2: coefficeint of determination

% Get time
et0 = clock;

% 0- check basic data dimension requirements
[n,p] = size(X);
if n < 3
    disp('Error: PRESS_GLM: data dimension is lower than 3!!')
    R2_cv = [];
    return
else
    if verbose
        fprintf(1,'PRESS_GLM: leave-one-out cross-validation scheme:...')
        fprintf(1,'%6.2f %%',0)
    end
end

% 1- perform leave-one-out scheme
Eg = NaN(size(y));
ii = 0;
for i=1:n
    X0 = X;
    X0(i,:) = [];
    y0 = y;
    y0(i,:) = [];
    b = pinv(X0'*X0)*X0'*y0;
    Eg(i,:) = X(i,:)*b;
    if verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',floor(100*i/n))
    end
end
if verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
    fprintf(1,'\n')
end

% 2- evaluate PRESS and R2_cv
for i=1:size(y,2)
    PRESS = sum((y(:,i)-Eg(:,i)).^2);
    TSS(i) = sum((y(:,i)-mean(y(:,i))).^2);
    R2_cv(i) = 1 - PRESS./TSS(i);
end

% 3- fit the model conditional on all data and wrap-up
if verbose
    disp('PRESS_GLM: model fit on full-data...')
end
iX = pinv(X'*X);
mP = iX*X'*y;
for i=1:size(y,2)
    yhat = X*mP(:,i);
    RSS = sum((yhat-y(:,i)).^2);
    vhat = RSS./(n-p);
    C = vhat*iX;
    VP(:,i) = diag(C);
    R2(i) = 1 - RSS./TSS(i);
end
if verbose
    disp('PRESS_GLM: model fit on full-data... OK.')
end
