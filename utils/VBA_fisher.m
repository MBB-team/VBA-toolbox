function [z] = VBA_fisher (r, inverse)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [z] = VBA_fisher (r, inv)
% apply Fisher transformation to Pearson correlation coefficient
%
% IN:
%   - r: correlation to transform
%   - inverse: flag to apply the reverse transformation (default = false)
%
% OUT:
%   - C: nxn correlation matrix
%
% /////////////////////////////////////////////////////////////////////////

if nargin < 2 || ~ inverse
    assert(VBA_isInRange(r, [-1 1]), ...
      '*** VBA_fisher: correlation coefficient should be between -1 and 1'); 
  
    z = 0.5 * log ((1 + r) ./ (1 - r));
else
    z = (exp(2 * r) - 1 ) ./ (1 + exp(2 * r));
end