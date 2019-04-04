function A = VBA_orth (A, norm)
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

if VBA_isWeird(A) 
    error('***VBA_orth: your matrix should not include any weird values (Nans, Inf, etc).');
end

[n,p] = size(A);
I = eye(n);

if norm
    A(:,1) = VBA_zscore (A(:, 1));
end
for i=2:p
    X = A(:, 1 : i - 1);
    nX = I - X * pinv(X' * X) * X';
    A(:,i) = nX * A(:,i);
    if norm
        A(:,i) = VBA_zscore (A(:,i));
    end
end
