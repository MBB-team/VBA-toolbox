function [flag] = VBA_isBinary (X)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [flag] = VBA_isBinary (X)
% check if X only contains 0 or 1 values
%
% IN:
%   - X: matrix, cell array of matrices, or structure to be checked
% OUT:
%   - flag: - true if all elements of X are binary (0 or 1)
%           - false otherwise
%           - NaN if any element or field is not numeric
%
% /////////////////////////////////////////////////////////////////////////

switch class (X)
    case 'logical'
        flag = true;
        
    case 'double'
        flag = all (ismember (X(:), [0, 1]));   
        
    case 'cell'
        flag = all (cellfun (@VBA_isBinary, X));
        
    case 'struct'
        flag = all (structfun (@VBA_isBinary, X));
    
    otherwise
        flag = NaN;
end


