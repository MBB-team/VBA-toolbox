function [sY] = smooth2(Y,va)
% smooth 2D image with gaussian kernel
% function [sY] = smooth2(Y,va)
% IN:
%   - Y: the input image
%   - va: the smoothing value

if ~exist('va','var') || isempty(va)
    va = 0;
end

if va > 0
    nv = 4;
    w = zeros(2*nv+1);
    for i=-nv:nv
        for j=-nv:nv
            w(i+nv+1,j+nv+1) = exp(-0.5.*(i.^2+j.^2)/va);
        end
    end
    w = w./sum(w(:));
    sY = conv2(Y,w,'same');
else
    sY = Y;
end
