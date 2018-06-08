function [ha,hf,lr] = unwrapVBvolatileOTO(posterior,out)
% OTO: recovers posterior sufficient statistics from VB inversion

if isempty(posterior)
    ha =[];
    return
end

[n,n_t] = size(posterior.muX);
if n==5
    
    ka = out.options.inF.lev2* VBA_sigmoid(posterior.muTheta(1), 'scale', out.options.inF.kaub);
    om = posterior.muTheta(2);
    mux = posterior.muX([2,4],:);
    sx = exp(posterior.muX([3,5],:));
    
    hf = figure('color',[1 1 1]);
    ha = subplot(2,2,1,'parent',hf,'ygrid','on','xlim',[1,n_t]);
    plotUncertainTimeSeries(mux,sx,1:n_t,ha);
    title(ha,'Learner''s belief estimate')
    legend(ha,{'x2 = outcome probability [LOG space]','x3 = expected volatility [LOG space]'})
    box(ha,'off')
    
    ha(2) = subplot(2,2,2,'parent',hf,'nextplot','add','ygrid','on','xlim',[1,n_t]);
    plot(ha(2),VBA_sigmoid(mux(1,:)));
    mv = exp(ka*mux(2,:)+om);
    plot(ha(2),mv,'g');
    
    % plot belief about outcome probability
    er1 = mux(1,:) - sqrt(sx(1,:));
    er2 = mux(1,:) + sqrt(sx(1,:));
    yp = [VBA_sigmoid(er2),fliplr(VBA_sigmoid(er1))];
    xp = [1:size(mux,2),fliplr(1:size(mux,2))];
    hfi = fill(xp,yp,'r','facecolor','b','edgealpha',0,'facealpha',0.25,'parent',ha(2));
    ev1 = exp(ka*(mux(2,:)-sqrt(sx(2,:)))+om);
    ev2 = exp(ka*(mux(2,:)+sqrt(sx(2,:)))+om);
    yp = [ev2,fliplr(ev1)];
    hfi = fill(xp,yp,'r','facecolor','g','edgealpha',0,'facealpha',0.25,'parent',ha(2));
    legend(ha(2),{'outcome probability','volatility'})
    box(ha(2),'off')
    
    
    x = posterior.muX;
    x(3,:) = exp(x(3,:));
    x(5,:) = exp(x(5,:));
    s1h = VBA_sigmoid(x(2,:)).*(1-VBA_sigmoid(x(2,:))); % likelihood precision
    s2h = x(3,:) + exp(ka*x(4,:)+om); % 2nd-level prediction variance
    lr = 1./(s2h.^-1 + s1h); % posterior variance
    ha(3) = subplot(2,2,3,'parent',hf,'nextplot','add','ygrid','on','xlim',[1,n_t]);
    plot(ha(3),lr);
    title(ha(3),'LEARNING RATE')
    box(ha(3),'off')
    
else
    
    ka = out.options.inF.lev2*VBA_sigmoid(posterior.muTheta(1),'scale',out.options.inF.kaub);
    om = posterior.muTheta(2);
    mux = posterior.muX([2,4,7,9],:);
    sx = exp(posterior.muX([3,5,8,10],:));
    
    hf = figure('color',[1 1 1]);
    ha = subplot(2,2,1,'parent',hf,'ygrid','on','xlim',[1,n_t]);
    plotUncertainTimeSeries(mux,sx,1:n_t,ha);
    title(ha,'Learner''s belief estimate (log-space)')
    legend(ha,{'1st cue: outcome probability','1st cue: expected volatility','2nd cue: outcome probability','2nd cue: expected volatility'})
    box(ha,'off')
    
    ha(2) = subplot(2,2,2,'parent',hf,'nextplot','add','ygrid','on','xlim',[1,n_t]);
    mr = mux([1,3],:);
    mv = exp(ka*mux([2,4],:)+om);
    col = getColors(4);
    for i=1:2
        plot(ha(2),VBA_sigmoid(mr(i,:)),'color',col(2*(i-1)+1,:));
        plot(ha(2),mv(i,:),'color',col(2*(i-1)+2,:));
    end
    er1 = mux([1,3],:) - sqrt(sx([1,3],:));
    er2 = mux([1,3],:) + sqrt(sx([1,3],:));
    ev1 = exp(ka*(mux(2,:)-sqrt(sx(2,:)))+om);
    ev2 = exp(ka*(mux(2,:)+sqrt(sx(2,:)))+om);
    for i=1:2
        yp = [VBA_sigmoid(er2(i,:)),fliplr(VBA_sigmoid(er1(i,:)))];
        xp = [1:size(mux,2),fliplr(1:size(mux,2))];
        hfi = fill(xp,yp,'r','facecolor',col(2*(i-1)+1,:),'edgealpha',0,'facealpha',0.25,'parent',ha(2));
        yp = [ev2,fliplr(ev1)];
        hfi = fill(xp,yp,'r','facecolor',col(2*(i-1)+2,:),'edgealpha',0,'facealpha',0.25,'parent',ha(2));
    end
    legend(ha(2),{'1st cue: outcome probability','1st cue: expected volatility','2nd cue: outcome probability','2nd cue: expected volatility'})
    box(ha(2),'off')
    title(ha(2),'Learner''s belief estimate')
    
    x = posterior.muX;
    x([3,5,8,10],:) = exp(x([3,5,8,10],:));
    s1h = VBA_sigmoid(x([2,7],:)).*(1-VBA_sigmoid(x([2,7],:))); % likelihood precision
    s2h = x([3,8],:) + exp(ka*x([4,9],:)+om); % 2nd-level prediction variance
    lr = 1./(s2h.^-1 + s1h); % posterior variance
    ha(3) = subplot(2,2,3,'parent',hf,'nextplot','add','ygrid','on','xlim',[1,n_t]);
    plot(ha(3),lr');
    title(ha(3),'LEARNING RATES')
    box(ha(3),'off')
    
end


X = [posterior.muX];
C = corrcoef(X');
ha(4) = subplot(2,2,4,'parent',hf);
imagesc(C,'parent',ha(4))
axis(ha(4),'square')
colorbar('peer',ha(4))
set(ha(4),'clim',[-1,1])
title(ha(4),'STATES'' CORRELATION')

VBA_getSubplots ();


