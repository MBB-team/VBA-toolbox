function [posterior,out] = VBA_PRESS(y,u,f_fname,g_fname,dim,options)

% Train/test predictive power evaluation for VBA model inversion.
% function [posterior,out] = VBA_PRESS(y,u,f_fname,g_fname,dim,options)
% This function evaluates the predictive power of any nonlinear state-space
% model of the form:
%   y_t = g( x_t,u_t,phi ) + e_t
%   x_t = f( x_t-1,u_t,theta ) + f_t
% using a VBA scheme (see VBA_NLStateSpaceModel.m).
% The deafult train/test cross-validation procedure is a leave-one-out
% scheme, and the predictive power is measured in terms of the resulting
% percentage of variance explained (R2_cv). First, one computes the
% predictive sum-of-square residual error (PRESS) as follows:
% PRESS = sum_i (y_i-E[g(x)|y_\i])^2
% where E[g(x)|y_\i] is the model prediction for y_i, given all the data
% but y_i.
% The predictive power is then defined as:
% R2_cv = 1 - PRESS/TSS
% where TSS = sum_i (y_i-E[y])^2 is the total data variance.
% IN: [see VBA_NLStateSpaceModel.m]. Additional inputs are provided
% through the structure options.cv, which contains the following fields:
%   .verbose: verbose mode
%   .DisplayWin: 
% OUT: [see VBA_NLStateSpaceModel.m]. Additional outputs are provided
% through the structure out.cv, which contains the following fields:
%   .press: the predictive sum-of-square residual error
%   .R2: the predictive power, in terms of R2_cv



try
    options.isYout;
catch
    options.isYout = zeros(size(y));
end

try
    options.cv.verbose;
catch
    options.cv.verbose = 1;
end

try
    options.verbose;
catch
    options.verbose = 0;
end

try
    options.DisplayWin;
catch
    options.DisplayWin = 0;
end


% Get time
et0 = clock;

% 0- check basic data dimension requirements
indIn = find(options.isYout==0);
if length(indIn) < 3
    disp('Error: VBA_PRESS: data dimension is lower than 3!!')
    posterior = [];
    out = [];
    return
else
    if options.cv.verbose
        fprintf(1,'VBA_PRESS: leave-one-out cross-validation scheme:...')
        fprintf(1,'%6.2f %%',0)
    end
end

% 1- perform leave-one-out scheme
Eg = NaN(size(y));
ii = 0;
for i=1:size(y,1)
    for j=1:size(y,2)
        options0 = options;
        if options.isYout(i,j)==0
            options0.isYout(i,j) = 1;
            [p,o] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options0);
            Eg(i,j) = o.suffStat.gx(i,j);
            ii = ii +1;
            if options.cv.verbose
                fprintf(1,repmat('\b',1,8))
                fprintf(1,'%6.2f %%',floor(100*ii/length(indIn)))
            end
        end
    end
end

% Display progress
if options.verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
    fprintf(1,'\n')
end

% 2- evaluate PRESS and R2_cv
PRESS = sum((y(indIn)-Eg(indIn)).^2);
TSS = sum((y(indIn)-mean(y(indIn))).^2);
R2_cv = 1 - PRESS./TSS;

% 3- invert the model conditional on all data and wrap-up
VBA_disp({' ','VBA_PRESS: model inversion given full-data...'},options.cv)
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
out.cv.press = PRESS;
out.cv.R2 = R2_cv;
VBA_disp('VBA_PRESS: model inversion given full-data... OK.',options.cv)

