function [] = displaySimulations(y,x,eta,e)
% plots simulated time series (including state-space SVD projections)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
if isweird({y,x,eta,e})
    return
end

dTime = [1:size(x,2)];

hf = figure;
pos = get(hf,'position');
set(hf,...
    'name','Simulated time series',...
    'position',[pos(1),pos(2)-pos(4),pos(3),2*pos(4)],...
    'color',ones(1,3),...
    'menubar','none');


ha = subplot(3,2,1,'parent',hf);
plot(ha,dTime,x')
grid(ha,'on')
axis(ha,'tight')
xlabel(ha,'time')
ylabel(ha,'x')
title('simulated hidden-states time series')
ha = subplot(3,2,2,'parent',hf,'nextplot','add');
plot(ha,dTime,y',':')
plot(ha,dTime,y','.')
plot(ha,dTime,y'-e')
% plot(ha,dTime,y'-e','.')
grid(ha,'on')
axis(ha,'tight')
xlabel(ha,'time')
ylabel(ha,'y')
title('simulated observations')
if size(x,1) > 3
    [u,s,v] = svd(x,0);
    xp = u(1:3,1:3)*s(1:3,:)*v';
elseif size(x,1) == 1
    xp = [x;zeros(2,size(x,2))];
elseif size(x,1) == 2
    xp = [x;zeros(1,size(x,2))];
else
    xp = x;
end
ha = subplot(3,2,3,'parent',hf);
plot3(ha,xp(1,:),xp(2,:),xp(3,:),'.')
set(ha,'nextplot','add')
plot3(ha,xp(1,:),xp(2,:),xp(3,:))
grid(ha,'on')
axis(ha,'tight')
title(ha,'x: state space')
xlabel(ha,'x1')
ylabel(ha,'x2')
zlabel(ha,'x3')


if size(y,1) > 3
    [u,s,v] = svd(y,0);
    yp = u(1:3,1:3)*s(1:3,:)*v';
elseif size(y,1) == 1
    yp = [y;zeros(2,size(y,2))];
elseif size(y,1) == 2
    yp = [y;zeros(1,size(y,2))];
else
    yp = y;
end
ha = subplot(3,2,4,'parent',hf);
plot3(ha,yp(1,:),yp(2,:),yp(3,:),'.')
set(ha,'nextplot','add')
plot3(ha,yp(1,:),yp(2,:),yp(3,:))
grid(ha,'on')
axis(ha,'tight')
title(ha,'y: state space')
xlabel(ha,'y1')
ylabel(ha,'y2')
zlabel(ha,'y3')

ha = subplot(3,2,5);
plot(ha,dTime,eta);
grid(ha,'on')
axis(ha,'tight')
title(ha,'stochastic innovations')
xlabel(ha,'time')
ylabel(ha,'eta')

ha = subplot(3,2,6);
plot(ha,dTime,e);
grid(ha,'on')
axis(ha,'tight')
title(ha,'measurement noise')
xlabel(ha,'time')
ylabel(ha,'e')
getSubplots

