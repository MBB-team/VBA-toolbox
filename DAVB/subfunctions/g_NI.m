function [gx] = g_NI(x,phi,u,in)
% dummy non-identifiable observation function
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
try
    X = in.X;
catch
    try
        X = ones(in.n,1);
    catch
        X = ones(10,1);
    end
end

gx = X*(phi(1)+phi(2));
