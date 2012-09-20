function handles = VBA_displayGroupBMC(posterior,out)

[K,n] = size(out.L);

hf = findobj('tag','groupBMC');
if isfield(out.options,'handles') && ismember(out.options.handles.hf,hf)
    
    handles = out.options.handles;
    
else
    
    pos0 = get(0,'screenSize');
    pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.85*pos0(4)];
    handles.hf = figure(...
        'position',pos,...
        'color',[1 1 1],...
        'name','group-level Bayesian model comparison',...
        'tag','groupBMC');
    handles.ha(1) = subplot(3,2,1,'parent',handles.hf,'nextplot','add');
    handles.ha(2) = subplot(3,2,2,'parent',handles.hf,'nextplot','add');
    handles.ha(3) = subplot(3,2,3,'parent',handles.hf,'nextplot','add');
    handles.ha(4) = subplot(3,2,4,'parent',handles.hf,'nextplot','add');
    handles.ha(5) = subplot(3,2,5,'parent',handles.hf,'nextplot','add');

    % display data (ie log-model evidences)
    col = getColors(n);
    for i=1:n
        plot(handles.ha(1),out.L(:,i),'color',col(i,:))
        plot(handles.ha(1),out.L(:,i),'.','color',col(i,:))
    end
    xlabel(handles.ha(1),'models')
    ylabel(handles.ha(1),'log- model evidences')
    set(handles.ha(1),...
        'xtick',1:K,...
        'xlim',[0.5,K+0.5],...
        'ygrid','on')
    title(handles.ha(1),'log- model evidences')
    
end



% display model attributions
cla(handles.ha(2))
hi = imagesc(posterior.r','parent',handles.ha(2));
axis(handles.ha(2),'tight')
% axis(handles.ha(2),'equal')
xlabel(handles.ha(2),'models')
ylabel(handles.ha(2),'subjects')
title(handles.ha(2),'model attributions')
set(handles.ha(2),...
    'xlim',[0.5,K+0.5],...
    'xtick',[1:K],...
    'ylim',[0.5,n+0.5],...
    'ytick',[1:n],...
    'clim',[0 1])
colorbar('peer',handles.ha(2))


% display model frequencies
cla(handles.ha(3))
[haf,hf,hp] = plotUncertainTimeSeries(out.Ef,diag(out.Vf),[],handles.ha(3));
% bar(handles.ha(3),out.Ef,'facecolor',0.5*[1 1 1])
plot(handles.ha(3),[0.5,K+0.5],[1,1]/K,'r')
xlabel(handles.ha(3),'models')
ylabel(handles.ha(3),'estimated model frequencies')
set(handles.ha(3),...
    'xtick',1:K,...
    'xlim',[0.5,K+0.5],...
    'ylim',[0 1],...
    'ygrid','on')
title(handles.ha(3),'estimated model frequencies')

% display model frequencies
cla(handles.ha(4))
bar(handles.ha(4),out.ep,'facecolor',0.5*[1 1 1])
plot(handles.ha(4),[0.5,K+0.5],[0.95,0.95],'r')
xlabel(handles.ha(4),'models')
ylabel(handles.ha(4),'exceedance probabilities')
set(handles.ha(4),...
    'xtick',1:K,...
    'xlim',[0.5,K+0.5],...
    'ylim',[0 1],...
    'ygrid','on')
title(handles.ha(4),'approximated exceedance probabilities')

% display VB algorithm convergence
cla(handles.ha(5))
plot(handles.ha(5),out.F,'k')
plot(handles.ha(5),out.F,'k.')
[haf,hf,hp] = plotUncertainTimeSeries(out.F0.*[1,1],[3,3].^2,[0.5,length(out.F)+0.5],handles.ha(5));
set(hp,'color',[1 0 0])
set(hf,'facecolor',[1 0 0])
% plot(handles.ha(5),[0.5,length(out.F)+0.5],out.F0.*[1,1],'r')
xlabel(handles.ha(5),'VB iterations')
ylabel(handles.ha(5),'VB free energy')
set(handles.ha(5),...
    'xtick',1:length(out.F),...
    'xticklabel',[],...
    'xlim',[0.5,length(out.F)+0.5],...
    'ygrid','on')
title(handles.ha(5),'VB algorithm convergence')

try;getSubplots;end

drawnow


