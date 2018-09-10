function [a, b] = VBA_guessHyperpriors (y, rangeEV)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [a, b] = VBA_guessHyperpriors (y, rangeEV)
% propose shape and scale parameters to use as hyperpriors for the
% measurement noise precision, based on the part of the variance one can
% expect the model to explain.
%
% IN:
%   - y: data to be fitted
%   - rangeEV: range of the fraction of variance to be explained. 
%       default = [0.1, 0.9], .ie between 10% and 90% of the total variance
% OUT:
%   - a, b: shape and scale parameters of the Gamma hyperprior
%
% /////////////////////////////////////////////////////////////////////////


%% check parameters
% =========================================================================
% fill in defaults
if nargin == 1
    rangeEV = [0.1, 0.9];
end

% check inputs
assert (numel (y) > 3, '*** VBA_guessHyperpriors: no enough datapoints!');
assert (numel (rangeEV) == 2, '*** VBA_guessHyperpriors: rangeEV must have 2 values.');
assert (VBA_isInRange (rangeEV, [1e-8, 1 - 1e-8]), '*** VBA_guessHyperpriors: rangeEV must be between 0 and 1.');

% cleanup
y(isnan (y)) = [];

%% Compute expected precision
% =========================================================================
% reformulate statistics of interest
varianceTotal = var (VBA_vec (y));
varianceUnexplained = varianceTotal * (1 - rangeEV);
precisionExpected = 1 ./ varianceUnexplained;

% shortcuts
sumPrecision = sum (precisionExpected);
diffPrecision = max (precisionExpected) - min (precisionExpected);

% precision range as 98% confidence interval
a = 6 * (sumPrecision / diffPrecision) ^ 2;
b = 12 * sumPrecision / (diffPrecision ^ 2);

