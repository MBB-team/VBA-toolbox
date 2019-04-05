function [y, out] = VBA_removeOutliers (y, t)

% // VBA toolbox //////////////////////////////////////////////////////////
%
% [y, out] = VBA_removeOutliers (y, t)
% remove outliers based on robust z-scoring
%
% IN:
%   - y: data to diagnostic
%   - t: threshold for the outlier detection (default = 3.5) 
%
% OUT:
%   - y: data cleaned of outliers
%   - out: flags indicating which points in original y were positively
%       detected as outliers
%
% This function computes robust estimates of the two first moments of the
% sample's distribution, using the median and the median absolute
% deviation (MAD). Outliers are defined as the samples that deviate more
% than t standard deviations from the first moment.
%
% By default, the function apply conservative threshold of 3.5 following
% Iglewicz, B., & Hoaglin, D. C. (1993). How to detect and handle outliers
%
% /////////////////////////////////////////////////////////////////////////

% check input parameters
% ----------------------------------------------
if nargin < 2
    t = 3.5;
end

% median absolute deviation
% ----------------------------------------------
my = median (y);
MAD = median (abs (y - my));

% consistent estimator of the standard deviation
% ----------------------------------------------
% scale factor for the normal distribtion, 
% b = 1 / Q(0.75), P.J. Huber (1981) Robust statistics
b = 1.4826; 
% estimator
sigmaHat = b * MAD;

% standardize data
% ----------------------------------------------
z = abs (y - my) / sigmaHat;

% threshold
% ----------------------------------------------
out = z > t;
y(out) = [];
