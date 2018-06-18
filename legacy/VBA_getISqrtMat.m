function S = VBA_getISqrtMat(C,inv)
% legacy code
s = warning ('on');
warning ('*** The function `VBA_getISqrtMat` is now deprecated. Please see `VBA_sqrtm` for an alternative.') 
warning (s);

% fallback
if nargin < 2
    inv = true;
end
S = VBA_getISqrtMat(C, inv);