function K = mykron(A,B)
%KRON Kronecker tensor product.
% KRON(X,Y) is the Kronecker tensor product of X and Y.
% The result is a large matrix formed by taking all possible
% products between the elements of X and those of Y. For
% example, if X is 2 by 3, then KRON(X,Y) is
%
% [ X(1,1)*Y X(1,2)*Y X(1,3)*Y
% X(2,1)*Y X(2,2)*Y X(2,3)*Y ]
%
% If either X or Y is sparse, only nonzero elements are multiplied
% in the computation, and the result is sparse.

% Jordan Rosenthal, 12/02/99, revision to original Matlab KRON function:
% Paul L. Fackler, North Carolina State, 9-23-96
% Copyright (c) 1984-98 by The MathWorks, Inc.
% $Revision: 5.10 $ $Date: 1997/11/21 23:39:55 $

[ma,na] = size(A);
[mb,nb] = size(B);



   ia = 1:ma;
   ia = ia(ones(mb,1),:);
   ib = (1:mb)';
   ib = ib(:,ones(ma,1));
   ja = (1:na);
   ja = ja(ones(nb,1),:);
   jb = (1:nb)';
   jb = jb(:,ones(na,1));
   K = A(ia,ja).*B(ib,jb);

