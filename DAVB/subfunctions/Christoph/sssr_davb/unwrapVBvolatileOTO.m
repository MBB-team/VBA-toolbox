function [ha] = unwrapVBvolatileOTO(posterior,out)
% OTO: recovers posterior sufficient statistics from VB inversion

if isempty(posterior)
    ha =[];
    return
end

[n,n_t] = size(posterior.muX);

mux = posterior.muX([2,4],:);
sx = posterior.muX([3,5],:);

hf = figure('color',[1 1 1]);
ha = subplot(3,1,1);
set(ha,'ygrid','on')
plotUncertainTimeSeries(mux,sx,1:n_t,ha);
title(ha,'Learner''s belief estimate')
legend(ha,{'x2','x3'})
box(ha,'off')

ha(2) = subplot(3,1,2);
set(ha(2),'nextplot','add','ygrid','on')
er1 = mux(1,:) - sqrt(sx(1,:));
er2 = mux(1,:) + sqrt(sx(1,:));
plot(ha(2),sigm(mux(1,:)));
yp = [sigm(er2),fliplr(sigm(er1))];
xp = [1:size(mux,2),fliplr(1:size(mux,2))];
hf = fill(...
    xp,yp,'r',...
    'facecolor','b',...
    'edgealpha',0,...
    'facealpha',0.25,...
    'parent',ha(2));
title(ha(2),'sgm( x2 ) = p(x1=1)')
box(ha(2),'off')

ka = posterior.muTheta(1);
om = posterior.muTheta(2);
ha(3) = subplot(3,1,3);
set(ha(3),'nextplot','add','ygrid','on')
er1 = mux(2,:) - sqrt(sx(2,:));
er2 = mux(2,:) + sqrt(sx(2,:));
plot(ha(3),exp(ka*mux(2,:)+om));
yp = [exp(ka*er2+om),fliplr(exp(ka*er1+om))];
xp = [1:size(mux,2),fliplr(1:size(mux,2))];
hf = fill(...
    xp,yp,'r',...
    'facecolor','b',...
    'edgealpha',0,...
    'facealpha',0.25,...
    'parent',ha(3));
title(ha(3),'exp( ka*x3 + om ) = STD[x2]')
box(ha(3),'off')

getSubplots