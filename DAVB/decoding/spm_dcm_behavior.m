function [] = spm_dcm_behavior(DCM)
    
    if ~isstruct(DCM)
        error('*** P should be a DCM structure');
    end
    
    %== compute profile
    fprintf('Preprocessing...\n');
    DCM.sensitivityProfile=SensitivityProfile(DCM.VBA.posterior,DCM.VBA.out);

    pos0 = get(0,'screenSize');
    pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];

    hfp = figure(...
    'position',pos,...
    'color',[1 1 1],...
    'name','DCM - Behavior',...
    'Renderer','OpenGL',...
    'menu','none');

    [ud.tabs.handles] = spm_uitab(hfp,{'Inputs',},{@contrast});
set(ud.tabs.handles.htab,'backgroundcolor',[1 1 1])
set(ud.tabs.handles.hh,'backgroundcolor',[1 1 1])
set(ud.tabs.handles.hp,'HighlightColor',0.8*[1 1 1])
set(ud.tabs.handles.hp,'backgroundcolor',[1 1 1])

ud.DCM = DCM;
set(hfp,'userdata',ud);

feval(@contrast,hfp)


end




function contrast(hfp)

ud=get(hfp,'userdata');

stims=ud.DCM.sensitivityProfile.stims;

all=horzcat(ud.DCM.sensitivityProfile.stims{:});
ulim_min = min(all,[],2);
ulim_max = max(all,[],2);
du = .05*(ulim_max-ulim_min);

n_s=length(stims);
n_u = size(stims{1},1);

for s=1:n_s
    stim=stims{s};
    for u=1:n_u
        li = sub2ind([n_u,n_s], u, s);
    subplot(n_s,n_u,li,...
        'parent',ud.tabs.handles.hp,...
        'visible','off');
    plot(stim(u,:)')
    set(gca,'box','off');
    set(gca,'XTick',[]);
    set(gca,'YTick',[ulim_min(u) ulim_max(u)]);
    set(gca,'YLim',[ulim_min(u)-du(u) ulim_max(u)+du(u)]);
    if s==1;
        title(sprintf('u_%d',u));
    end

    set(gca,'visible','on');
    end

end

end