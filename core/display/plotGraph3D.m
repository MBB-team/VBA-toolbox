function [hi,px,gx,x] = plotGraph3D(px,gx,x,ha,va,norm)
% plot image as a map in 3D
% function [hi] = plotGraph3D(px,gx,x,ha,va)
% IN:
%   - px: the pxn image to plot
%   - gx: the nX1 grid on the first dimension of px
%   - x: a 1Xp trajectory to plot onto the 3D surface
%   - ha: the handle of the parent axes
%   - va: the smoothing value

if ~exist('x','var')
    x = [];
end
if ~exist('ha','var') || isempty(ha)
    hf = figure('color',[1 1 1]);
    ha = axes('parent',hf);
end
if ~exist('va','var') || isempty(va)
    va = 0;
end
if ~exist('norm','var') || isempty(norm)
    norm = 1;
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
    px = conv2(px,w);
    dx = mean(diff(gx));
    offset = dx*[1:nv];
else
    offset = [];
    nv = 0;
end

if norm % normalize px
    px =  bsxfun(@rdivide,px,sum(px,2));
end
mig = min(gx);
mag = max(gx);
gx = [mig-flipud(offset(:));gx;mag+offset(:)];
gt = [-nv:size(px,1)-nv-1]';
hi = surf(gt,gx,px','parent',ha);
set(hi,'edgecolor',0.8*[1 1 1],'edgealpha',0.1),

if isempty(x) % compute first-order moment from time-dependent density
   x = zeros(1,length(gt)-2*nv);
   for t=1:length(gt)-2*nv
       pxt = px(t+nv,:);
       x(t) = sum(pxt.*gx');
   end
end

% add deterministic trajectory
[ny,nx] = hist(repmat(x,2,1),gx);
ny = [zeros(length(gx),nv),ny,zeros(length(gx),nv)];
ny = double(~~ny);
ny(ny==0) = NaN;
z = px'.*ny;
zz = z(~isnan(z));
xx = x(~isnan(x));
tt = repmat(gt',size(px,2),1).*ny;
tt = tt(~isnan(tt));
set(ha,'nextplot','add')
hh = plot3(ha,tt,xx,zz,'linewidth',2);

set(ha,'nextplot','replace')
axis(ha,'tight')
box(ha,'on')
xlabel('time')
ylabel('value')
zlabel('probability density')
colormap(flipud(autumn))
try, VBA_getSubplots (); end

