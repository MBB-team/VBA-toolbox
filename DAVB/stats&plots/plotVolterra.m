function [ha] = plotVolterra(hf,mw,vw)

nu = size(mw,2);

if nu==6
    
    str={'choice=1 * feedback',...
        'choice=0 * feedback',...
        'opponent=1 * feedback',...
        'opponent=0 * feedback',...
        'choice=1 * [CST]',...
        'opponent=1 * [CST]'};
    
    n = nu/2+1;
    for i=1:nu
        ha(i) = subplot(n,2,i,'parent',hf);
        plotUncertainTimeSeries(mw(:,i),vw(:,i),[],ha(i));
        title(ha(i),str{i})
        xlabel(ha(i),'lag')
        ylabel(ha(i),'weight')
    end
    
    ha(nu+1) = subplot(n,2,2*n-1,'parent',hf);
    plotUncertainTimeSeries(mw(:,1)-mw(:,2),vw(:,1)+vw(:,2),[],ha(nu+1));
    xlabel(ha(nu+1),'lag')
    ylabel(ha(nu+1),'weight')
    title(ha(nu+1),'[net effect of feedback]')
    
    ha(nu+2) = subplot(n,2,2*n,'parent',hf);
    plotUncertainTimeSeries(mw(:,3)-mw(:,4),vw(:,3)+vw(:,4),[],ha(nu+2));
    xlabel(ha(nu+2),'lag')
    ylabel(ha(nu+2),'weight')
    title(ha(nu+2),'[net effect of opponent''s feedback]')
    
elseif nu==4
    
    str={'own feedback',...
        'opponent''s feedback',...
        'own action',...
        'opponent''s action'};
    
    n = nu/2;
    for i=1:nu
        ha(i) = subplot(n,2,i,'parent',hf);
        plotUncertainTimeSeries(mw(:,i),vw(:,i),[],ha(i));
        title(ha(i),str{i})
        xlabel(ha(i),'lag')
        ylabel(ha(i),'weight')
    end
    
elseif nu==3
    
    str={'own feedback',...
        'opponent''s feedback',...
        'own action'};
    
    n = 2;
    for i=1:nu
        ha(i) = subplot(n,2,i,'parent',hf);
        plotUncertainTimeSeries(mw(:,i),vw(:,i),[],ha(i));
        title(ha(i),str{i})
        xlabel(ha(i),'lag')
        ylabel(ha(i),'weight')
    end
    
elseif nu==2
    
    str={'own action',...
        'opponent''s action'};
    
    n = nu/2;
    for i=1:nu
        ha(i) = subplot(n,2,i,'parent',hf);
        plotUncertainTimeSeries(mw(:,i),vw(:,i),[],ha(i));
        title(ha(i),str{i})
        xlabel(ha(i),'lag')
        ylabel(ha(i),'weight')
    end
    
    
end

