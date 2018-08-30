function [val, varargout] = VBA_maxMat (A)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [maxVal, idx1, idx2, ..., idxn] = VBA_maxMat (A)
% Find maximal value and its position in the matrix A. 
%
% IN:
%   - A: n-dimensional matrix

% OUT:
%   - val: maximal value in A
%   - idx1, idx2, ...: Position of the maximum value in A. 
%       * If only one index is requested ([val, idx] = VBA_maxMat (A))
%         then the linear position is returned (ie. A(idx) = val).
%       * Otherwise, the n indices specify the position of val in the
%         respective dimensions of A
%
% If the maximum value occurs more than once, then VBA_maxMat returns the
% index on the first occurence only
%
% /////////////////////////////////////////////////////////////////////////

% get maximum value and linear index
% -------------------------------------------------------------------------
[val, idx] = max (A(:));

% find position in matrix
% -------------------------------------------------------------------------
s = size(A);
switch nargout
    % no position
    case 1
        return
    % linear index
    case 2
        varargout{1} = idx;
    % subscript indices 
    case numel (s) + 1
        varargout = cell (size (s));
        [varargout{:}] = ind2sub (s, idx);
    % wrong number of output!
    otherwise
        error ('*** VBA_maxMat the number of output indices should match the input dimension.')
end