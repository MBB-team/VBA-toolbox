function variance = nanvar2(y)
% Calculate the variance of a vector by removing nan-values.
% IN:
%   - y: The data, can be a matrix or a vector.
% OUT:
%   - variance: The variance of y. 

y = y(isnan(y) == 0);
variance = var(y);

end