function iA = pinvComplex(A)
% computes complex matrix inversion
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
if any(imag(A(:)))
    
    rA = real(A);
    cA = imag(A);
    n = size(A,1);
    AA = [  rA      cA
            -cA     rA ];
    iAA = VB_inv(AA);
    iA = iAA(1:n,1:n) + sqrt(-1).*iAA(1:n,n+1:2*n);
else
    iA = VB_inv(A);
end
