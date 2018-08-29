function [handles] = VBA_projectMoG(posterior,out,y)
% display clustering results on 2D data eigenspace
% function [handles] = plotResults(y,xi,eta,F,theta,K_opt,options)
% This function projects the MoG components onto the first two
% eigencomponents of the data. This 2D eigenspace is then tesselated
% according to which component is the most likely at each point in space.
% This allows for a rough (visual) determination of the components'
% boundaries, in the 2D data eigenspace.
% IN:
%   - posterior: structure containing the following fields:
%       .muEta: pxK prior mean of the components 1st-order moment
%       .SigmaEta: Kx1 cell array of posterior pxp covariance matrices of
%       the components 1st-order moment
%       .a_gamma: Kx1 shape parameter of the posterior on the components'
%       2nd-order moment
%       .a_gamma: Kx1 scale parameter of the posterior on the components'
%       2nd-order moment
%       .d: Kx1 vector of posterior counts for each component
%       .Z: nxK matrix containing the estimated probability, for each
%       voxel, of belonging to each of the clusters
%   - out:
%       .dt: elapsed time (in sec)
%       .dim: final dimension of the model
%       .options: (filled in) options for MoG VB model inversion
%       .suffStat: internal sufficient statistics
%       .date: date (vector format)
%       .F: history of Free Energy across VB iterations
%       .it: # VB iterations until convergence
%   - y: the pxn data matrix, where n is the number of 'observations' to be
%   classified, and p is the dimension of each of these observations
% OUT:
%   - handles: the handles of the output graphical objects

[p,n] = size(y);
K = size(posterior.muEta,2);
if p < 2
    handles = [];
    return
end

% 0- initialize display window
try
    set(out.handles.ha(9),'nextplot','add')
    handles.ha = out.handles.ha(9);
    set(handles.ha,'xlim',[1,1e2],'ylim',[1,1e2])
catch
    handles.hf = figure('name','Classification result: 2D eigen-projection','color',[1 1 1]);
    handles.ha = axes('parent',handles.hf,'nextplot','add');
    title(handles.ha,['Final number of components: K=',num2str(K),' , log p(y|K=',num2str(K),')>',num2str(out.F(end),'%4.3e')])
    xlabel(handles.ha,'first data eigen-mode (%)')
    ylabel(handles.ha,'second data eigen-mode (%)')
end
% 1- identify 2D data eigenspace
if p > 2
    [u,s,v,] = svd(y,0);
    y = s(1:2,:)*v';
else
    u = eye(2);
end
% 2- Evaluate the data likelihood of each MoG's component on a 100x100 grid 
miy = min(y,[],2);
may = max(y,[],2);
dy = (may-miy).*1e-2;
[X,Y] = meshgrid(miy(1):dy(1):may(1),miy(2):dy(2):may(2));
lpdf = zeros(length(X),length(Y),K);
mu = zeros(p,K);
muk = zeros(2,K);
for k = 1:K
    mu(:,k) = posterior.muEta(:,k);
    vark = posterior.b_gamma(k)./posterior.a_gamma(k);
    if out.options.normalize
        mu(:,k) = diag(sqrt(diag(out.normalize.Q)))*mu(:,k) + out.normalize.m;
        vark = vark.*out.normalize.Q;
    end
    muk(:,k) = u(:,1:2)'*mu(:,k);
    vark = u(:,1:2)'*vark*u(:,1:2);
    [lpdf(:,:,k)] = getLogNormpdf(muk(:,k),vark,X,Y);
end
% 3- Detect and display components' boundaries
[mp,imp] = max(lpdf,[],3);
[Ix,Iy] = gradient(imp);
di = abs(Ix) + abs(Iy);
[C,h] = contour(handles.ha,~~di,1);
set(h,'color',[0.5 0.5 0.5])
% 4- project and plot each data point on the 2D eigenspace
y(1,:) = y(1,:) - miy(1);
y(2,:) = y(2,:) - miy(2);
y(1,:) = y(1,:)./(dy(1));
y(2,:) = y(2,:)./(dy(2));
[me,im] = max(posterior.z',[],2);
col = getColors(K);
for i = 1:n
    plot(handles.ha,y(1,i),y(2,i),'.','color',col(im(i),:));
end
% 5- project and plot the components' modes
muk = (muk - repmat(miy(:),1,K))./repmat(dy(:),1,K);
for k = 1:K
    text(muk(1,k),muk(2,k),num2str(k),'color',col(k,:),'parent',handles.ha,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
end
grid(handles.ha,'on')
try;VBA_getSubplots ();end

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
