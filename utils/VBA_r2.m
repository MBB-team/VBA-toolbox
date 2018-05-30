function r2 = VBA_r2 (pred, data, isYout)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% r2 = VBA_r2 (pred, data, isYout)
% coefficient of determination
%
% IN:
%   - pred: predictions of the classifier
%   - data: predicted data. 
%   - isYout: array of flags to exclude datapoints (default = no exclusion)
%
% OUT:
%   - r2: coefficient of determination (fraction of explained variance)
%
% /////////////////////////////////////////////////////////////////////////

% check input parameters
% =========================================================================

% default to no exclusions
if nargin < 3
    isYout = zeros (size (pred));
end

% ensure exclusion
pred(logical (isYout)) = nan;
data(logical (isYout)) = nan;

% coefficient of determination
% =========================================================================

SS2_tot = sumall ((data(:) - VBA_nanmean(data(:))) .^2);
SS2_err = sumall ((data(:) - pred(:)) .^2);

r2 = 1 - (SS2_err / SS2_tot);

end

function s = sumall(z)
    s = sum (sum (z, 'omitnan'), 'omitnan');
end