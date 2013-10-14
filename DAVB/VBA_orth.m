function A = VBA_orth(A,norm)
% serial orthogonalisation of matrix A
% function A = VBA_orth(A,norm)
% This function considers each column serially, projecting it in the null
% space of the previous columns of A.
% IN:
%   - A: nXp matrix
%   - norm: a flag for normalization of A's columns
% OUT:
%   - A: the serially orthogonalized matrix A
try,norm;catch,norm=0;end
[n,p] = size(A);
I = eye(n);
if norm && std(A(:,1))~=0
    A(:,1) = (A(:,1)-mean(A(:,1)))./std(A(:,1));
end
for i=2:p
    X = A(:,1:i-1);
    nX = I - X*pinv(X'*X)*X';
    A(:,i) = nX*A(:,i);
    if norm && std(A(:,i))~=0
        A(:,i) = (A(:,i)-mean(A(:,i)))./std(A(:,i));
    end
end
