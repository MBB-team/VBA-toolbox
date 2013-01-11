function [handles] = plotResults(y,xi,eta,F,theta,K_opt,options)

if ~exist('options','var') || isempty(options) || ~isfield(options,'verbose')
    options.verbose = 1;
end
[p,n] = size(y);
if p < 2 || ~options.verbose
    handles = [];
    return
end

handles.hf = figure('name','Classification result: 2D eigen-projection','color',[1 1 1]);
handles.ha = gca(handles.hf);
title(handles.ha,['Final number of components: K=',num2str(K_opt),' , log p(y|K=',num2str(K_opt),')>',num2str(F(end),'%4.3e')])
xlabel(handles.ha,'first data eigen-mode (%)')
ylabel(handles.ha,'second data eigen-mode (%)')
set(handles.ha,'nextplot','add')
if p > 2 % project onto first two eigenmodes
    [u,s,v,] = svd(y,0);
    y = s(1:2,:)*v';
else
    u = eye(2);
end
miy = min(y,[],2);
may = max(y,[],2);
dy = (may-miy).*5e-3;
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
LPDF = reshape(lpdf,size(lpdf,1)*size(lpdf,2),[]);
[mp,imp] = max(LPDF,[],2);
I = zeros(size(lpdf,1),size(lpdf,2));
for k=1:K_opt
    I(imp==k) = k;
end
[Ix,Iy] = gradient(I);
di = abs(Ix) + abs(Iy);
[C,h] = contour(handles.ha,~~di,1);
set(h,'color',[0.5 0.5 0.5])
y(1,:) = y(1,:) - miy(1);
y(2,:) = y(2,:) - miy(2);
y(1,:) = y(1,:)./(dy(1));
y(2,:) = y(2,:)./(dy(2));
% get most probable label for each observation
[me,im] = max(xi,[],2);
col = getColors(K_opt);
for i = 1:n
    plot(handles.ha,y(1,i),y(2,i),'.','color',col(im(i),:));
end
muk = u(:,1:2)'*eta;
muk = (muk - repmat(miy(:),1,K_opt))./repmat(dy(:),1,K_opt);
for k = 1:K_opt
    plot(handles.ha,muk(1,k),muk(2,k),'+','color',col(k,:));
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
