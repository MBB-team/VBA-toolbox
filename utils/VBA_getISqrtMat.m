function S = VBA_getISqrtMat(C,inv)
% quick derivation of the (inverse) square root of a matrix
% function S = VBA_getISqrtMat(C,inv)
% IN:
%   - C: the entry matrix
%   - inv: binary flag for inverse square root (inv=1) or square root
%   (inv=0) matrix operator
% OUT:
%   - S: the (inverse) square root of C.

if nargin < 2
    inv = 1;
else
    inv = ~~inv;
end
C(C==Inf) = realmax;  % dirty fix for infinite precision matrices
if sum(C(:)) ~= 0
    if isequal(C,diag(diag(C)))
        if inv
            S = diag(sqrt(diag(C).^-1));
        else
            S = diag(sqrt(diag(C)));
        end
    else
        [U,s,V] = svd(full(C));
        if inv
            S = U*diag(sqrt(diag(s).^-1))*V';
        else
            S = U*diag(sqrt(diag(s)))*V';
        end
    end
else
    S = 0;
end
