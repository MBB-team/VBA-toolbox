function [z, mu, sigma] = VBA_zscore (x, w, dim, nanflag)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [z, mu, sigma] = VBA_zscore (x, flag, dim, nanflag)
% standard score, aka z-score
%
% Compute the stardardized zscore of x, ie z = (x - mean(x)) / std(x)
%
% /////////////////////////////////////////////////////////////////////////

if nargin < 2 || isempty (w)
	w = 0;
end

if nargin < 3 || isempty (dim)
	dim = find (size(x) > 1, 1);
	if isempty(dim)
		dim = 1;
	end
end

if nargin < 4
    nanflag = 'includenan';
end

mu = mean (x, dim, nanflag);
sigma = std (x, w, dim, nanflag);
sigma(sigma == 0) = 1;

z = bsxfun(@minus, x, mu);
z = bsxfun(@rdivide, z, sigma);
