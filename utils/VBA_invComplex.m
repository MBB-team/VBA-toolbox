function iA = VBA_invComplex (A)

% // VBA toolbox //////////////////////////////////////////////////////////
%
% iA = VBA_invComplex (A)
% complex matrix inversion
%
% IN:
%   - A: matrix to invert, potentially complex
%
% OUT:
%   - iA: inverse of A
%
% /////////////////////////////////////////////////////////////////////////

if isreal (A)
    iA = VBA_inv(A);
    return
end

rA = real (A);
cA = imag (A);
n = size (A, 1);
    
AA = [  rA, cA
      - cA, rA ];
iAA = VBA_inv (AA);
iA = iAA(1 : n, 1 : n) + sqrt(- 1) .* iAA(1 : n, n + 1 : 2 * n);

