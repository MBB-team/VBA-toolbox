function [accuracy, balanced_accuracy] = VBA_accuracy (pred, data, type, isYout)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [accuracy, balanced_accuracy] = VBA_accuracy (pred, data, type, isYout)
% classification scores of binary or categorical predictor 
%
% IN:
%   - pred: predictions of the classifier
%   - data: predicted data. 
%   - type: flag for binary (1) or categorical (2) data (default = 1)
%   - isYout: array of flags to exclude datapoints (default = no exclusion)
%
% Note that for categorical data (k classes), predictions are kx1 vectors
% (columns of pred) summing to one and data are kx1 indicator vectors (all
% 0 but for the observed class)
%
% OUT:
%   - accuracy: proportion of correct classification 
%   - balanced_accuracy: accuracy corrected for chance level
%
% Note that balanced accuracyis only computed for binary data (so far)
%
% /////////////////////////////////////////////////////////////////////////

% check input parameters
% =========================================================================
% default to binomial
if nargin < 3
    type = 1;
end

% catch gaussian case
if type == 0 
    accuracy = NaN;
    balanced_accuracy = NaN;
    return
end

% default to no exclusions
if nargin < 4
    isYout = zeros (size (pred));
end

% ensure exclusion
pred(logical (isYout)) = nan;
data(logical (isYout)) = nan;

% find best guess
% =========================================================================
% number of classes
switch type
    case 1 % binomial
        k = 2;       
    case 2 % multinomial
        k = size (pred, 1);
    otherwise
        error('*** VBA_accuracy: type must be 1 or 2');
end

% binarized model predictions
bg = pred > (1 / k);  

% simple accuracy
% =========================================================================
% ratio of best guesses matching the data

% true positives
tp = sumall (data .* bg);
% false positives
fp = sumall ((1 - data) .* bg);
% false positives
fn = sumall (data .* (1 - bg));
% true negatives
tn = sumall ((1 - data) .* (1 - bg));

% total positive predictions
P = tp + fn;
% total negative predictions
N = tn + fp;

% classification score
accuracy = (tp + tn) ./ (P + N);

% balanced accuracy
% =========================================================================
% correct score for chance level

switch type
    case 1 % binomial
        balanced_accuracy = 0.5 * (tp ./ P + tn ./ N);
        
    case 2 % multinomial
        balanced_accuracy = NaN;
        % TODO: find a good score!
        % Arzheimer & Evans (2013) A New Multinomial Accuracy Measure for Polling Bias, Political Analysis
end

end

function s = sumall(z)
    s = nansum (nansum (z));
end