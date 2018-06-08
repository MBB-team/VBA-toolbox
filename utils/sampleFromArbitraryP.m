function [X] = sampleFromArbitraryP(p,gridX,N)
% inverse transform sampling scheme
% function [X] = sampleFromArbitraryP(p,gridX,N)
% This function samples from an arbitrary 1D probability distribution
% IN:
%   - p: pX1 vector (the density evaluated along the grid)
%   - gridX: pX1 vector (the grid over which the density is evaluated)
%   - N: the number of samples to be sampled
% OUT:
%   - X: NX1 vector of samples

try; N; catch, N=1; end
p = VBA_vec(p);

if size(gridX,1)==1
    gridX = VBA_vec(gridX);
end
k = size(gridX,2);

if any(isnan(p))
    X = nan(N,k);
else
    pcdf = cumsum(p(:));
    X = zeros(N,k);
    for i=1:N
        below = find(rand<=pcdf);
        X(i,:) = gridX(below(1),:);
    end
end