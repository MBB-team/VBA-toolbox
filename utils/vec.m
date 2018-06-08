function vx = vec(X)
% computes the Vec operator
% function vx = vec(X)
% JD, 2/03/2007.

if isempty(X)
    vx = [];
else
    vx = full(X(:));
end


