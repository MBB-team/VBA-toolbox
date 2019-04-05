function [z, mu, sigma] = VBA_zscore (x, flag, dim, omitnan)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [z, mu, sigma] = VBA_zscore (x, flag, dim, ignore_nans)
% standard score, aka z-score
%
% Compute the stardardized zscore of x, ie z = (x - mean(x)) / std(x)
%
% IN:
%   x: vector or matrix to normalize
%   flag: if true, normalizes x using std(X,1), i.e., by computing the
%     standard deviation(s) using N rather than N-1, where N is the length of
%     the dimension along which VBA_zscore works. VBA_zscore(x,0) is the same as
%     VBA_zscore(x). Default = false;
%   dim: dimension along which the standardization is performed. By
%     default, uses the first non singleton dimension.
%   ignore_nans: if true, use nan-safe normalization (default = false).
%
% OUT: z-score of x
%
% /////////////////////////////////////////////////////////////////////////

if nargin < 2 || isempty (flag)
	flag = 0;
else
    flag = +(flag > 0);
end

if nargin < 3 || isempty (dim)
	dim = find (size(x) > 1, 1);
	if isempty(dim)
		dim = 1;
	end
end

if nargin < 4
    omitnan = false;
end

if omitnan
    mu = VBA_nanmean (x, dim);
    sigma = sqrt(VBA_nanvar (x, dim, flag));
else
    mu = mean (x, dim);
    sigma = std (x, flag, dim);
end
sigma(sigma == 0) = 1;

z = bsxfun(@minus, x, mu);
z = bsxfun(@rdivide, z, sigma);
