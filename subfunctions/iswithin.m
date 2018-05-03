function [flag] = iswithin(X,bounds)
% true if all elements of X is within bounds
% function [flag] = iswithin(X, bounds)
% IN:
%   - X: N-D matrix (or cell array of matrices) to be checked
%   - bounds: a 2X1 vector that defines an interval on the real line
% OUT:
%   - flag: 1 if all elements of X are within the bounds, 0 if not

if VBA_isWeird (X)
    flag = 0;
    return
end

flag = 1;
if iscell(X)
    for i=1:numel(X)
        flag = flag & iswithin(X{i},bounds);
    end
elseif isstruct(X)
    fn = fieldnames(X);
    for i=1:length(fn)
        flag = flag & iswithin(getfield(X,fn{i}),bounds);
    end
elseif isnumeric(X) || islogical(X)
    if any(X(:)<min(bounds)) || any(X(:)>max(bounds))
        flag = 0;
    end
end