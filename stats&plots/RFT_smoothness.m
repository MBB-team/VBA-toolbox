function [fwhm] = RFT_smoothness(X)
% estimates the smoothness of a Gaussian RF sampled on a regular lattice
% [fwhm] = RFT_smoothness(X)
% Note: if X is a Lxn matrix, then RFT_smoothness assumes that columns of X
% are repeated measures of the field (whose dimension is L).
% IN:
%   - X: Lxn matrix of sampled RF
% OUT:
%   - fwhm: field's smoothness (in voxels)

[L,n] = size(X);
if n>1 % normalize field across repetitions (if #repetitions>1)
    sx2= var(X,[],2);
    X = X./repmat(sqrt(sx2),1,n);
end
dX = (X(1:L-2,:) - X(3:L,:))./2;
fwhm = sqrt(4*log(2)./var(dX(:)));