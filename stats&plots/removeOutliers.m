function [yc,out] = removeOutliers(y)
% removes outlies based upon robust moment-matched Gaussian distribution
% function [yc,out] = removeOutliers(y)
% The logic of outlier removal is as follows. Let z~N(0,1) be a standard
% Gaussian variable. If we sample n times the parent distribution N(0,1),
% then the expected number of events [z>t], where t is some arbitrary
% threshold would then be given by n*P(z>t). If the expected number of such
% events is below one, then we do not expect any value in our sample to
% exceed t. In other word, any value that exceeds t (such that n*P(z>t) <1)
% is an outlier.
% This function computes robust estimates of the two first moments of the
% sample's parent distribution, using the median and the median absolute
% deviation (MAD). Outliers are defined as those samples, for which the
% expected number of more extreme events is below unity.
% IN:
%   - y: nx1 data sample
% OUT:
%   - yc: (n-k)x1 corrected sample, having removed k outliers
%   - out: the indices of the detected outliers

n = length(y);
my = median(y);
sy = 1.4826*median(abs(y-my));
z = abs(y-my)./sy;
En = n*(1-spm_Ncdf(z,0,1));
out = find(En<1);
yc = y;
yc(out) = [];
