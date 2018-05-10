function [dfdx, fx] = numericDiff(fName, idxArg2Diff, varargin)
% legacy code
s = warning ('on');
warning ('*** The function `numericDiff` is now deprecated. Please use `VBA_numericDiff` (same syntax).') 
warning (s);

% fallback
[dfdx, fx] = VBA_numericDiff(fName, idxArg2Diff, varargin)