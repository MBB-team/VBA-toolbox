function y = VBA_nanmean(x, dim)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% y = VBA_nanmean(varargin)
% mean ignoring NaNs
%
% This function enhances the functionality of NANMEAN as distributed in
% the MATLAB Statistics Toolbox and is meant as a replacement (hence the
% identical name).  
%
% NANMEAN(X,DIM) calculates the mean along any dimension of the N-D
% array X ignoring NaNs.  If DIM is omitted NANMEAN averages along the
% first non-singleton dimension of X.
%
% See also MEAN
% -------------------------------------------------------------------------
%    author:      Jan Glaescher
%    affiliation: Neuroimage Nord, University of Hamburg, Germany
%    email:       glaescher@uke.uni-hamburg.de
%    
%    $Revision: 1.1 $ $Date: 2004/07/15 22:42:13 $
%
% /////////////////////////////////////////////////////////////////////////

if isempty(x)
	y = NaN;
	return
end
if nargin < 2
	dim = min(find(size(x)~=1));
	if isempty(dim)
		dim = 1;
	end
end
% Replace NaNs with zeros.
nans = isnan(x);
x(isnan(x)) = 0; 
% denominator
count = size(x,dim) - sum(nans,dim);
% Protect against a  all NaNs in one dimension
i = find(count==0);
count(i) = ones(size(i));
y = sum(x,dim)./count;
y(i) = i + NaN;
% $Id: nanmean.m,v 1.1 2004/07/15 22:42:13 glaescher Exp glaescher $
