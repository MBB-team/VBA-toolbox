function S = VBA_weighted_inner_product (A, B, w)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% S = VBA_weighted_inner_product (A, B, w)
% computes inner products between 
%
% IN:
%   - A: matrix or column vector of real values
%   - B: matrix or column vector of real values (optional, if not specified
%        or empty, it is set to A)
%   - w: column vector of real values (optional, if not specified or empty
%        it is set to 1)
%        or matrix of real values
%   
% OUT:
%   - S: Weighted inner product between column vectors in A and B given the
%        weight vector w, such that S(i,j) = Sum_k A(k,i) w(k) B(k,j) or in
%        matlab notation sum(A(:,i) .* w .*  B(:, j)).
%        If w is a matrix (supposed to be symmetric positive definite, but
%        this is not checked), then S = A' * W * B.
%
% /////////////////////////////////////////////////////////////////////////

if nargin<2 || isempty(B)
    B = A;
end
if nargin<3 || isempty(w)
    w = 1;
end


assert(size(A,1)==size(B,1), ...
    'VBA:invalidInput', ...
    '*** VBA_weighted_inner_product: A and B must have the same inner dimensions');


if isvector(w)
    assert(length(w)==1 || length(w)==size(A,1),  ...
        'VBA:invalidInput', ...
        '*** VBA_weighted_inner_product: Weight vector w must have length 1 or the same dimension as A and B');
    
    S = A' * bsxfun(@times, w(:), B);
else
    assert(ismatrix(w) && all(size(w)==size(A,1)),  ...
        'VBA:invalidInput', ...
        '*** VBA_weighted_inner_product: Weight matrix w must have dimension N x N with N being the first dimension of A and B');
    
    S = A' * w * B;
end
