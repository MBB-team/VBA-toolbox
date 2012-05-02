function [Q,ha] = risk_OTO(posterior,out)

% OTO: inferred speed-accuracy risk from reaction times

alpha = exp(posterior.muPhi(1));
theta = exp(posterior.muPhi(2));

% dlambda0 = posterior.muX(end,:);

mat = max(out.y);
mit = max([0,min(out.y)]);
dt = (mat-mit)*1e-3;
gridt = mit:dt:mat;
gridl0 = -0.999:1e-3:1;
dyn = 1-exp(-2.*alpha*gridt);

Q = zeros(length(gridl0),length(gridt));
opt = zeros(length(gridl0),2);
c = 1; % posterior risk is symmetric wrt to c and sign(dlambda0)!
for i=1:length(gridl0)
    lt = gridl0(i).*dyn;
    Q(i,:) = lt.*(1-2*c) + c + theta.*gridt;
    [tmp,miQ] = min(Q(i,:));
    opt(i,:) = [gridl0(i),miQ];
end

hf = figure('color',ones(1,3));
pos = get(hf,'position');
set(hf,'position',pos+[0,-pos(4),0,pos(4)])
ha = subplot(2,1,1,'parent',hf);
set(ha,'nextplot','add')
imagesc(Q,'parent',ha)
axis(ha,'tight')
colorbar('peer',ha,'location','NorthOutside');

xl = get(ha,'xlim');
yl = get(ha,'ylim');
% reduce to 0-1:
opt(:,1) = (opt(:,1)-min(opt(:,1)))./((max(opt(:,1))-min(opt(:,1))));
opt(:,2) = (opt(:,2)-min(opt(:,2)))./((max(opt(:,2))-min(opt(:,2))));
% rescale to [xl,yl]
opt(:,1) = diff(yl).*opt(:,1) + yl(1);
opt(:,2) = diff(xl).*opt(:,2) + xl(1);
contour(ha,Q,'color',0.5*ones(1,3))
plot(ha,xl,mean(yl)*ones(1,2),...
    'color',0.5*ones(1,3),...
    'linestyle','--')
plot(ha,...
    opt(:,2),opt(:,1),...
    'color',ones(1,3),...
    'linewidth',2)
ytick = get(ha,'ytick');
for i=1:length(ytick)
    ytl{i} = num2str(gridl0(ytick(i)));
end

set(ha,...
    'ytick',ytick,...
    'yticklabel',ytl)
xlabel(ha,'time (sec)')
ylabel(ha,'posterior-prior expectation')
title(ha,'posterior risk (subject''s choice = 1)')
box(ha,'on')

if any(out.y<0)
    disp('Warning: negative reation times!')
    out.y(out.y<0) = 0;
end
[ny,nx] = hist(out.y);
ha(2) = subplot(2,1,2,'parent',hf);
hb = bar(nx,1e2*ny./sum(ny),'parent',ha(2),'facecolor',[.8 .8 .8]);
set(ha(2),'xlim',[min(out.y),mat]);
xlabel(ha(2),'time (sec)')
ylabel(ha(2),'frequency (%)')
title(ha(2),'empirical distribution of reaction time data')
grid(ha(2),'on')

xtick = get(ha(2),'xtick');
xtick = (xtick-min(out.y))./((mat-min(out.y)));
xtick = diff(xl).*xtick + xl(1);
xtl = get(ha(2),'xticklabel');
set(ha(1),...
    'xtick',xtick,...
    'xticklabel',xtl)





