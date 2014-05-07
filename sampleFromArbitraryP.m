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
p = vec(p);
if size(gridX,1)==1
    gridX = vec(gridX);
end
pcdf = cumsum(p(:));
k = size(gridX,2);
X = zeros(N,k);
for i=1:N
    below = find(rand<=pcdf);
    X(i,:) = gridX(below(1),:);
end