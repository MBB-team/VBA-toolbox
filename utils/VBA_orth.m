function A = VBA_orth (A, norm, nanflag)
% serial orthogonalisation of matrix A
% function A = VBA_orth(A,norm)
% This function considers each column serially, projecting it in the null
% space of the previous columns of A.
% IN:
%   - A: nXp matrix
%   - norm: a flag for normalization of A's columns
% OUT:
%   - A: the serially orthogonalized matrix A

if nargin < 2
    norm = false;
end

if nargin < 3
    nanflag = 'includenan';
end

[n,p] = size(A);
I = eye(n);

if norm
    A(:,1) = VBA_zscore (A(:, 1), [], [], nanflag);
end
for i=2:p
    X = A(:, 1 : i - 1);
    nX = I - X * pinv(X' * X) * X';
    A(:,i) = nX * A(:,i);
    if norm
        A(:,i) = VBA_zscore (A(:,i), [], [], nanflag);
    end
end
