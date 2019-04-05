function vx = VBA_vec(X)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% vx = VBA_vec (X)
% transform array X into a column vector
%
% IN:
%   - X: matrix
% OUT:
%   - vx: vector
%
% /////////////////////////////////////////////////////////////////////////

% avoid weirdly sized empty vector
if isempty(X)
    vx = [];
    return;
end

% return vectorized input
vx = full(X(:));


