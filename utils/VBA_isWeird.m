function [flag] = VBA_isWeird (X)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [flag] = VBA_isWeird (X)
% check if X contains any Infs, NaNs or non real entries
%
% IN:
%   - X: matrix, cell array of matrices, or structure to be checked
% OUT:
%   - flag: - true if any element or field of X is weird (Infs, Nans or non-real)
%           - false if no element or field contains 
%           - NaN if any element or field is not numeric
%
% /////////////////////////////////////////////////////////////////////////

switch class (X)
    case 'double'
        flag = any (isinf (X(:)) | isnan (X(:)) | ~ isreal (X(:)));   
        
    case 'cell'
        flag = any (cellfun (@VBA_isWeird, X));
        
    case 'struct'
        flag = any (structfun (@VBA_isWeird, X));
    
    otherwise
        flag = NaN;
end


