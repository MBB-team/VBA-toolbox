function hf = unwrapKTOM(x,inG)
% display k-ToM's evolving beliefs over trials
% function unwrapKTOM(x,options)
% IN:
%   - x: k-ToM's hidden-states
%   - inG: the input structure to g_kToM.m
% OUT:
%   - hf: display figure handle

K = inG.lev; % sophistication level

nt = size(x,2); % nb of trials
Ntot = inG.indParev + inG.indParobs; % total nb of params per opponent

hf = figure('color',[1 1 1],'name',[num2str(K),'-ToM learner']);

   
if K > 1
    
    % plot belief about opponent's sophistication
    na = 3;
    Pk = VBA_sigmoid(x(1:(K-1),:)); % P(k), with k=0,...,k'-2
    Pk = [Pk;1-sum(Pk,1)]; % insert last P(k=k'-1)
    ha(1) = subplot(na,1,1,'parent',hf,'nextplot','add');
    plot(ha(1),Pk')
    plot(ha(1),[1,nt],[1,1]./K,'color',0.5*[1 1 1])
    xlabel(ha(1),'trials')
    ylabel(ha(1),'P(op=k-ToM)')
    for k=1:K
        leg{k} = ['op = ',num2str(k-1),'-ToM'];
    end
    legend(ha(1),leg)
    title(ha(1),'opponent''s sophistication')
    set(ha(1),'xlim',[1,nt],'ylim',[0,1])
    
else
    na = 2;
end

% plot belief about opponent's next move
ha(2) = subplot(na,1,na-1,'parent',hf,'nextplot','add');
if K >0
    Po = NaN(K,nt);
    for k=1:K
        Po(k,:) = VBA_sigmoid(x(inG.indlev(k).f,1:nt));
    end
else
    Po = VBA_sigmoid(x(1,:));
end
plot(ha(2),Po')
plot(ha(2),[1,nt],[1,1]./2,'color',0.5*[1 1 1])
xlabel(ha(2),'trials')
if K >1
    ylabel(ha(2),'P(o=1|op=k-ToM)')
    legend(ha(2),leg)
elseif K==1
    ylabel(ha(2),'P(o=1|op=0-ToM)')
else
    ylabel(ha(2),'P(o=1)')
end
title(ha(2),'opponent''s next move')
set(ha(2),'xlim',[1,nt],'ylim',[0,1])



% plot belief about opponent's parameters
ha(na) = subplot(na,1,na,'parent',hf,'nextplot','add');
if K > 0
    mP = NaN(Ntot,nt,K); % E[params]
    vP = NaN(Ntot,nt,K); % V[params]
    for k=1:K
        indP = inG.indlev(k).Par;
        mP(:,:,k) = x(indP(1:2:2*Ntot),1:nt);
        vP(:,:,k) = exp(x(indP(2:2:2*Ntot),1:nt));
    end
else
    mP = x(1,:);
    vP = exp(x(2,:));
end
ud.K = K;
ud.mP = mP;
ud.vP = vP;
ud.ha = ha(na);
set(hf,'userdata',ud)
if K > 1
    hc = uicontrol('parent',hf,'style','popupmenu','string',leg,'callback',@DOdisplayOpParams);
    set(hc,'units','normalized','position',[0.75 0.03 0.2 0.05])
end
displayOpParams(ud,1)

try, VBA_getSubplots (); end

end



function DOdisplayOpParams(e1,e2)
op = get(e1,'value');
ud = get(get(e1,'parent'),'userdata');
displayOpParams(ud,op)
end

function displayOpParams(ud,op)
m = ud.mP(:,:,op);
v = ud.vP(:,:,op);
hc = get(ud.ha,'children');
delete(hc)
[haf,hf,hp] = plotUncertainTimeSeries(m,v,[1:size(m,2)],ud.ha);
xlabel(ud.ha,'trials')
if ud.K >1
    ylabel(ud.ha,'P(op params|op=k-ToM)')
elseif ud.K==1
    ylabel(ud.ha,'P(op params|op=0-ToM)')
else
    ylabel(ud.ha,'P(op params)')
end
title(ud.ha,'opponent''s parameters')
end
