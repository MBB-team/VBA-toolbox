function ldQ = VBA_logDet(Q,indIn,t)
% computes the log-determinant of matrix Q
% function ldQ = VBA_logDet(Q,indIn,t)
% IN:
%   - Q: nxn matrix
%   - indIn: vector of indices that defines a submatrix of Q, whose
%   determinant is computed {1:n}
%   - t: threshold on the eigenvalues of the submatrix {eps}
% OUT:
%   - ldQ: log-determinant of the submatrix of Q

if nargin < 3
    t = eps;
end
if nargin < 2 || isempty(indIn)
    dq = diag(Q);
    indIn = find(~isinf(dq)&dq~=0);
end
subQ = full(Q(indIn,indIn));
if isequal(subQ,eye(length(indIn))) % identity matrix
    ldQ = 0;
elseif isequal(subQ,diag(diag(subQ))) % diagonal matrix
    dQ = diag(subQ);
    ldQ = sum(log(dQ(abs(dQ)>t)));
else % full matrix
    if ~ VBA_isWeird (subQ)
        ev = eig(subQ);
    else
        ev = 0;
    end
    ev = ev(abs(ev)>t);
    ldQ = sum(log(ev));
end
ldQ = real(ldQ);

