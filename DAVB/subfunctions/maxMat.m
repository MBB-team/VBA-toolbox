function [iy,ix,Max] = maxMat(A)

% finds x- and y- indices of the maximum entry within a 2D-matrix
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

[my,iyx] = max(A,[],1);
[Max,ix] = max(my);
iy = iyx(ix);
