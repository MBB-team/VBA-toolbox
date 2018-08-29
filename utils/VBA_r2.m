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

% check inputs
narginchk(2,3);
assert (all (size (pred) == size (data)), 'VBA_r2: predictions and data must have the same size');
assert (all (size (isYout) == size (data)), 'VBA_r2: isYout and data must have the same size');

% catch trivial case
if isempty (data) || numel(data) < 3
    r2 = nan;
    return;
end

% ensure exclusion
pred(logical (isYout)) = nan;
data(logical (isYout)) = nan;

% coefficient of determination
% =========================================================================

SS2_tot = sumall ((data(:) - VBA_nanmean(data(:))) .^2);
SS2_err = sumall ((data(:) - pred(:)) .^2);

r2 = max (0, 1 - (SS2_err / SS2_tot));

end

function s = sumall(z)
    s = nansum (nansum (z));
end