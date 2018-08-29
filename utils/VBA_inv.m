function [iQ] = VBA_inv(Q,indIn,flag,v)
% overloaded sparse matrix pseudo-inverse
% function [iQ] = VBA_inv(Q,indIn,flag,v)
% IN:
%   - Q: the nxn matrix to be inverted
%   - indIn: a vector of indices that specifies the submatrix of Q that has
%   to be inverted. The rest of the matrix is padded with zeros. If empty,
%   the routine looks for infinite or below precision (close to zero)
%   entries in the the diagonal of Q.
%   - flag: if flag='replace', the routine returns Q, having replaced its
%   elements not in 'indIn' with v (see below)
%   - v: a number by which to padd the elements of Q not in 'indIn' (only
%   for flag='replace')
% OUT:
%   - iQ: the nxn matrix that is either the inverse of Q or v-padded Q (for
%   flag='replace').

if nargin < 2 || isempty(indIn)
    dq = diag(Q);
    indIn = find(~isinf(dq)&dq~=0);
end
% use lazy evaluation if all matrix is used
isSub = ~all(numel(indIn) == size(Q));

if nargin < 3
    replace = 0;
else
    replace = isequal(flag,'replace');
end
if nargin < 4
    v = 0;
end
if isSub
    subQ = full(Q(indIn,indIn));
else
    subQ = full(Q);
end
if replace % v-padd
    iQ = v.*ones(size(Q));
    iQ(indIn,indIn) = subQ;
else % (p)invert Q
    if isequal(subQ,eye(length(indIn)))   % identity matrix
        subiQ = subQ;
    elseif isequal(subQ,diag(diag(subQ))) % diagonal matrix
        tol  = max(eps(norm(diag(subQ),'inf'))*length(indIn),exp(-32)); 
        subiQ = diag((diag(subQ)+tol).^-1);
    else % full matrix
        tol  = max(eps(norm(subQ,'inf'))*length(indIn),exp(-32)); 
        subiQ = inv(subQ + eye(length(indIn))*tol);
    end
    if isSub
        iQ = zeros(size(Q));
        iQ(indIn,indIn) = subiQ;
    else
        iQ = subiQ;
    end
end



