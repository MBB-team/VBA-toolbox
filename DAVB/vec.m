function vx = vec(X)
% computes the Vec operator
% function vx = vec(X)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
% JD, 2/03/2007.

if isempty(X)
    vx = [];
else
    vx = full(X(:));
end


