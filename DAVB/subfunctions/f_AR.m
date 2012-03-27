function [fx,dfdx,dfdp,d2fdxdp] = f_AR(x,theta,u,in)
% AR(1) evolution function
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
fx = x;
dfdx = eye(length(x));
dfdp = [];
d2fdxdp = [];