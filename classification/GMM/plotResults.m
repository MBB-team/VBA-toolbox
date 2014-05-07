function [handles] = plotResults(y,xi,eta,F,theta,K_opt,options)
% display clustering results on 2D data eigenspace
% function [handles] = plotResults(y,xi,eta,F,theta,K_opt,options)
% This function projects the MoG components onto the first two
% eigencomponents of the data. This 2D eigenspace is then tesselated
% according to which component is the most likely at each point in space.
% This allows for a rough (visual) determination of the components'
% boundaries, in the 2D data eigenspace.
% IN:
%   - y : the pxn data matrix, where n is the number of
%   'observations' to be classified, and p is the dimension of each of
%   these observations
%   - xi: Kxn matrix containing the estimated probability, for each
%   observation, of belonging to each of the clusters
%   - eta: pxK prior mean of the components 1st-order moment
%   - F: tthe MoG log-evidence
%   - theta: structure containing the fields:
%       .gamma: the Kx1 vector of components' precisions
%   - K_opt: the number of components in the MoG
%   - options: structure containing the fields:
%       .verbose: flag for verbose mode (default=1)
%       .handles: grpahical objects handles
% OUT:
%   - handles: the handles of the output graphical objects

if ~exist('options','var') || isempty(options) || ~isfield(options,'verbose')
    options.verbose = 1;
end
    
if ~exist('options','var') || isempty(options) || ~isfield(options,'handles') || ~isfield(options.handles,'ha') || ~isequal(get(options.handles.ha,'type'),'axes')
    handles.hf = figure('name','Classification result: 2D eigen-projection','color',[1 1 1]);
    handles.ha = axes('parent',handles.hf);
    title(handles.ha,['Final number of components: K=',num2str(K_opt),' , log p(y|K=',num2str(K_opt),')>',num2str(F(end),'%4.3e')])
    xlabel(handles.ha,'first data eigen-mode (%)')
    ylabel(handles.ha,'second data eigen-mode (%)')
else
    handles = options.handles;
end

[p,n] = size(y);

set(handles.ha,'nextplot','add')

% 1- identify 2D data eigenspace
if p > 2
    [u,s,v,] = svd(y,0);
    y = s(1:2,:)*v';
elseif p == 2
    u = eye(2);
else
    handles = [];
    return
end
% 2- Evaluate the data likelihood of each MoG's component on a 100x100 grid 
miy = min(y,[],2);
may = max(y,[],2);
dy = (may-miy).*1e-2;
[X,Y] = meshgrid(miy(1):dy(1):may(1),miy(2):dy(2):may(2));
lpdf = zeros(length(X),length(Y),K_opt);
for k = 1:K_opt
    muk = eta(:,k);
    vark = theta.gamma(k).^-1;
    if p > 2
        muk = u(:,1:2)'*muk;
        vark = u(:,1:2)'*vark*u(:,1:2);
    end
    [lpdf(:,:,k)] = getLogNormpdf(muk,vark,X,Y);
end
% 3- Identify the most probable component at each voxel of the 100x100 grid
LPDF = reshape(lpdf,size(lpdf,1)*size(lpdf,2),[]);
[mp,imp] = max(LPDF,[],2);
I = zeros(size(lpdf,1),size(lpdf,2));
for k=1:K_opt
    I(imp==k) = k;
end
% 4- Detect and display components' boundaries
[Ix,Iy] = gradient(I);
di = abs(Ix) + abs(Iy);
[C,h] = contour(handles.ha,~~di,1);
set(h,'color',[0.5 0.5 0.5])
% 5- project and plot each data point on the 2D eigenspace
y(1,:) = y(1,:) - miy(1);
y(2,:) = y(2,:) - miy(2);
y(1,:) = y(1,:)./(dy(1));
y(2,:) = y(2,:)./(dy(2));
[me,im] = max(xi,[],2);
col = getColors(K_opt);
for i = 1:n
    plot(handles.ha,y(1,i),y(2,i),'.','color',col(im(i),:));
end
% 6- project and plot the components' modes
muk = u(:,1:2)'*eta;
muk = (muk - repmat(miy(:),1,K_opt))./repmat(dy(:),1,K_opt);
for k = 1:K_opt
    text(muk(1,k),muk(2,k),num2str(k),'color',col(k,:),'parent',handles.ha,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
end
grid(handles.ha,'on')


function lp = getLogNormpdf(muk,vark,X,Y)
ng = length(X);
lp = zeros(ng,length(Y));
[u,s,v] = svd(vark);
ds = diag(s);
lds = log(sum(ds(ds>0)));
ivark = pinv(vark);
for i=1:ng
    for j = 1:ng
        tmp = muk - [X(i,j);Y(i,j)];
        lp(i,j) = -0.5.*tmp'*ivark*tmp -0.5*lds - log(2*pi);
    end
end
