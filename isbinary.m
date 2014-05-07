function flag = isbinary(X)
% true if X contains only 0 or 1 entries
% function [flag] = isbinary(X)
% IN:
%   - X: N-D matrix (or cell array of matrices) to be checked
% OUT:
%   - flag: 1 id X is binary, 0 if not
if iscell(X)
    ok = 1;
    for i=1:numel(X)
        ok = ok & isbinary(X{i});
    end
    flag = ok;
elseif isstruct(X)
    ok = 1;
    fn = fieldnames(X);
    for i=1:length(fn)
        ok = ok & isbinary(getfield(X,fn{i}));
    end
    flag = ok;
elseif isnumeric(X) || islogical(X)
    flag = 0;
    if all(ismember(X(:),[0,1]))
        flag = 1;
    end
else
    flag = 0;
end

