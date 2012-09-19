function [gx] = g_NI(x,phi,u,in)

% dummy non-identifiable observation function

try
    X = in.X;
catch
    try
        X = ones(in.n,1);
    catch
        X = ones(10,1);
    end
end

gx = X*(sum(phi));
