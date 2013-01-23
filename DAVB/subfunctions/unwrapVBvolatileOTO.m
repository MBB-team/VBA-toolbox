function [ha] = unwrapVBvolatileOTO(posterior,out)
% OTO: recovers posterior sufficient statistics from VB inversion

if isempty(posterior)
    ha =[];
    return
end

[n,n_t] = size(posterior.muX);

mux = posterior.muX([2,4],:);
sx = posterior.muX([3,5],:);

ka = out.options.inF.lev2*sgm(posterior.muPhi(1),out.options.inF.kaub);
om = posterior.muPhi(2);
% th = sgm(posterior.muPhi(3),options.inF.thub);

hf = figure('color',[1 1 1]);
ha = subplot(3,1,1,'parent',hf,'ygrid','on','xlim',[1,n_t]);
plotUncertainTimeSeries(mux,sx,1:n_t,ha);
title(ha,'Learner''s belief estimate')
legend(ha,{'x2','x3'})
box(ha,'off')

ha(2) = subplot(3,1,2,'parent',hf,'nextplot','add','ygrid','on','xlim',[1,n_t]);
er1 = mux(1,:) - sqrt(sx(1,:));
er2 = mux(1,:) + sqrt(sx(1,:));
plot(ha(2),sigm(mux(1,:)));
yp = [sigm(er2),fliplr(sigm(er1))];
xp = [1:size(mux,2),fliplr(1:size(mux,2))];
hfi = fill(xp,yp,'r','facecolor','b','edgealpha',0,'facealpha',0.25,'parent',ha(2));
title(ha(2),'OUTCOME PROBABILITY: sgm( x2 ) = p(x1=1)')
box(ha(2),'off')

ha(3) = subplot(3,1,3,'parent',hf,'nextplot','add','ygrid','on','xlim',[1,n_t]);
er1 = mux(2,:) - sqrt(sx(2,:));
er2 = mux(2,:) + sqrt(sx(2,:));
plot(ha(3),exp(ka*mux(2,:)+om));
yp = [exp(ka*er2+om),fliplr(exp(ka*er1+om))];
xp = [1:size(mux,2),fliplr(1:size(mux,2))];
hf = fill(xp,yp,'r','facecolor','b','edgealpha',0,'facealpha',0.25,'parent',ha(3));
title(ha(3),'VOLATILITY: exp( ka*x3 + om ) = V[x2(t+1)|x(t)]')
box(ha(3),'off')

getSubplots