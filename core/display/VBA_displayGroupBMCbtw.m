function handles = VBA_displayGroupBMCbtw(ep,out)
% displays group-level BMC between-conditions analysis results
% function handles = VBA_displayGroupBMCbtw(ep,out)
% IN:
%   - ep/out: output structures of VBA_groupBMbtwC.m
% OUT:
%   - handles: structure containing the handles of the graphical objects




[nm,ns,nc] = size(out.L);

pos0 = get(0,'screenSize');
pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.85*pos0(4)];
handles.hf = figure('position',pos,'color',[1 1 1],'name','group-level Bayesian model comparison','tag','groupBMC');

% display data (ie log-model evidences)
handles.ha(1) = subplot(3,2,1,'parent',handles.hf,'nextplot','add');
col = getColors(nc);
leg = cell(nc,1);
for i=1:nc
    errorbar(handles.ha(1),mean(out.L(:,:,i),2),std(out.L(:,:,i),[],2)./sqrt(ns),'color',col(i,:),'marker','.')
    leg{i} = ['cond #',num2str(i)];
end
xlabel(handles.ha(1),'models')
set(handles.ha(1),'xtick',1:nm,'xlim',[0.5,nm+0.5],'ygrid','on')
VBA_title(handles.ha(1),'mean log- model evidences')
legend(handles.ha(1),leg)
drawnow

% display profile likelihood coherence
handles.ha(2) = subplot(3,2,2,'parent',handles.hf,'nextplot','add');
cla(handles.ha(2))
C = zeros(nc,nc,ns);
for n=1:ns
    for i=1:nc
        for j=1:nc
            tmp = corrcoef(VBA_vec(out.L(:,n,i)),VBA_vec(out.L(:,n,j)));
            C(i,j,n) = tmp(2,1);
        end
    end
end
hi = imagesc(mean(C,3),'parent',handles.ha(2));
axis(handles.ha(2),'tight')
axis(handles.ha(2),'square')
xlabel(handles.ha(2),'conditions')
ylabel(handles.ha(2),'conditions')
VBA_title(handles.ha(2),'profile likelihood coherence')
set(handles.ha(2),'xlim',[0.5,nc+0.5],'xtick',[1:nc],'ylim',[0.5,nc+0.5],'ytick',[1:nc],'ydir','reverse','clim',[-1 1])
colormap(handles.ha(2),bone)
handles.hc = colorbar('peer',handles.ha(2));
drawnow

if isfield(out.options,'families') && ~isempty(out.options.families)
    
    % display families partition
    nf = length(out.options.families);
    C = zeros(nm,nf);
    tmp = [];
    for i=1:nf
        indf = out.options.families{i};
        tmp = [tmp;VBA_vec(indf)];
        C(indf,i) = 1;
    end
    C = ~C; % for display purposes (black => model belongs to the family)
    handles.ha(3) = subplot(3,2,3,'parent',handles.hf,'nextplot','add');
    hi = imagesc(C','parent',handles.ha(3));
    axis(handles.ha(3),'tight')
    xlabel(handles.ha(3),'models')
    ylabel(handles.ha(3),'families')
    VBA_title(handles.ha(3),'families'' partition')
    set(handles.ha(3),'xlim',[0.5,nm+0.5],'xtick',[1:nm],'ylim',[0.5,nf+0.5],'ytick',[1:nf],'ydir','reverse','clim',[0 1])
    drawnow
    
    % display per-condition family-BMS
    C = zeros(nf,nc);
    for i=1:nc
        C(:,i) = VBA_vec(out.VBA.cond(i).out.families.ep);
    end
    handles.ha(5) = subplot(3,2,5,'parent',handles.hf,'nextplot','add');
    hi = imagesc(C','parent',handles.ha(5));
    axis(handles.ha(5),'tight')
    xlabel(handles.ha(5),'families')
    ylabel(handles.ha(5),'conditions')
    VBA_title(handles.ha(5),'per-condition family-BMS: EPs')
    set(handles.ha(5),'xlim',[0.5,nf+0.5],'xtick',[1:nf],'ylim',[0.5,nc+0.5],'ytick',[1:nc],'ydir','reverse','clim',[0 1])
    handles.hc(end+1) = colorbar('peer',handles.ha(5));
    drawnow
    
else
    
    % display per-condition BMS
    C = zeros(nm,nc);
    for i=1:nc
        C(:,i) = VBA_vec(out.VBA.cond(i).out.ep);
    end
    handles.ha(5) = subplot(3,2,5,'parent',handles.hf,'nextplot','add');
    hi = imagesc(C','parent',handles.ha(5));
    axis(handles.ha(5),'tight')
    xlabel(handles.ha(5),'models')
    ylabel(handles.ha(5),'conditions')
    VBA_title(handles.ha(5),'per-condition BMS: EPs')
    set(handles.ha(5),'xlim',[0.5,nm+0.5],'xtick',[1:nm],'ylim',[0.5,nc+0.5],'ytick',[1:nc],'ydir','reverse','clim',[0 1])
    handles.hc(end+1) = colorbar('peer',handles.ha(5));
    drawnow
    
end



if isfield(out,'factors')
    
    % display design
    sf = size(out.factors);
    sf(sf<=1) = [];
    nf = size(sf,2); % number of factors
    if nf > 1
        C = zeros(nc,nf);
        strpar = '[';
        for i=1:nf
            strpar = [strpar,'i',num2str(i),','];
        end
        strpar(end) = ']';
        for i=1:nc
            eval([strpar,'=ind2sub(size(out.factors),find(out.factors==i));'])
            for j=1:nf
                eval(['C(i,j)=i',num2str(j),';'])
            end
        end
        handles.ha(4) = subplot(3,2,4,'parent',handles.hf,'nextplot','add');
        hi = imagesc(C','parent',handles.ha(4));
        axis(handles.ha(4),'tight')
        xlabel(handles.ha(4),'conditions')
        ylabel(handles.ha(4),'factors')
        VBA_title(handles.ha(4),'design')
        set(handles.ha(4),'xlim',[0.5,nc+0.5],'xtick',[1:nc],'ylim',[0.5,nf+0.5],'ytick',[1:nf])
        handles.hc(end+1) = colorbar('peer',handles.ha(4));
        set(handles.hc(end),'ytick',1:max(VBA_vec(C)))
        drawnow
        
        % display pEP of per-factor btw-condition stability 
        pep = zeros(nf,1);
        for f=1:nf
            pep(f) = ep(f).*(1-out.VBA.btw(f).out.bor) + 0.5*out.VBA.btw(f).out.bor;
        end
        handles.ha(6) = subplot(3,2,6,'parent',handles.hf,'nextplot','add');
        hi = bar(ep,'parent',handles.ha(6));
        set(hi,'facecolor',0.8*ones(1,3))
        plot(handles.ha(6),pep,'ro')
        axis(handles.ha(6),'tight')
        xlabel(handles.ha(6),'factors')
        ylabel(handles.ha(6),'pEP(f=)')
        VBA_title(handles.ha(6),'btw-condition stability: pEPs')
        set(handles.ha(6),'xlim',[0.5,nf+0.5],'xtick',[1:nf],'ylim',[0,1])
        drawnow
        
    end
    
end

if ~isfield(out,'factors') || nf==1
    
    pep = ep.*(1-out.VBA.btw.out.bor) + 0.5*out.VBA.btw.out.bor;
    handles.ha(6) = subplot(3,2,6,'parent',handles.hf,'nextplot','add');
    hi = bar(ep,'parent',handles.ha(6));
    set(hi,'facecolor',0.8*ones(1,3))
    plot(handles.ha(6),pep,'ro')
    axis(handles.ha(6),'tight')
    xlabel(handles.ha(6),'f=')
    ylabel(handles.ha(6),'pEP(f=)')
    VBA_title(handles.ha(6),'btw-condition stability: pEPs')
    set(handles.ha(6),'xlim',[0.5,nf+0.5],'xtick',[],'ylim',[0,1])
    drawnow
    
end

try;VBA_getSubplots ();end




