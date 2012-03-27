function C = dMatdvec(A)

% Derivative of a matrix wrt to its non-zero elements
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

A = ~~A;
ind = find(A~=0);
n = numel(A);
ni = length(ind);
C = zeros(n,ni);
for i=1:length(ind)
    C(ind(i),i) = 1;
end