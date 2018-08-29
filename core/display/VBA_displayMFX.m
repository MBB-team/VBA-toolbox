function [options] = VBA_displayMFX(p_sub,o_sub,p_group,o_group,init,within)
% displays VBA_MFX inversion results

options = o_group.options;
try;init;catch;init=0;end

if ~options.DisplayWin
    return
end

try;within;catch;within='on';end


if init
    pos0 = get(0,'screenSize');
    pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];
    display.hfp = figure(...
        'position',pos,...
        'color',[1 1 1],...
        'menubar','none',...
        'name','VBA: MFX analysis',...
        'tag','VBA_MFX',...
        'Renderer','OpenGL',...
        'visible','off');
    for i=1:4
        for j=1:2
            display.ha(i,j) = subplot(4,2,(i-1)*2+j,'parent',display.hfp,'visible','off');
        end
    end
    display.ho = uicontrol(...
        'parent',display.hfp,...
        'style','text',...
        'units','normalized',...
        'position',[0.2,0.01,0.6,0.02],...
        'backgroundcolor',[1,1,1]);
    pos = get(display.ha(4,2),'position');
    str = cell(options.dim.ns,1);
    for i=1:options.dim.ns
        str{i} = ['subject #',num2str(i)];
    end
    display.ht = uicontrol(...
        'parent',display.hfp,...
        'style','text',...
        'string','MFX: within-subject VBA model inversions',...
        'units','normalized',...
        'position',[pos(1),pos(2)+0.12,pos(3),0.02],...
        'backgroundcolor',[1,1,1],...
        'visible',within);
    display.hl = uicontrol(...
        'parent',display.hfp,...
        'style','popupmenu',...
        'string',str,...
        'units','normalized',...
        'position',[pos(1),pos(2)+0.1,pos(3),0.02],...
        'backgroundcolor',[1,1,1],...
        'visible',within);
    display.hbg = uibuttongroup(...
        'parent',display.hfp,...
        'Position',[pos(1),pos(2)+0.07,pos(3),0.02],...
        'bordertype','none',...
        'backgroundcolor',[1,1,1],...
        'visible',within);
    display.hr1 = uicontrol(...
        'parent',display.hbg,...
        'style','radiobutton',...
        'string','final (empirical Bayes) priors',...
        'units','normalized',...
        'position',[0,0.1,0.6,0.8],...
        'backgroundcolor',[1,1,1],...
        'visible',within);
    display.hr2 = uicontrol(...
        'parent',display.hbg,...
        'style','radiobutton',...
        'string','initial priors',...
        'units','normalized',...
        'position',[0.7,0.1,0.3,0.8],...
        'backgroundcolor',[1,1,1],...
        'visible',within);
    display.hb = uicontrol(...
        'parent',display.hfp,...
        'style','pushbutton',...
        'string','DISPLAY selected subject',...
        'units','normalized',...
        'position',[pos(1),pos(2)+0.04,pos(3),0.02],...
        'backgroundcolor',[1,1,1],...
        'callback',@VBAMFX_displayWithin,...
        'visible',within);
    options.display = display;
    return
else
    try
        for i=1:4
            for j=1:2
                set(options.display.ha(i,j),'visible','off');
            end
        end
    catch
        [options] = VBA_displayMFX(p_sub,o_sub,p_group,o_group,1);
    end
    display = options.display;
    ud.p_sub = p_sub;
    ud.o_sub = o_sub;
    ud.o_group = o_group;
    ud.display = display;
    set(display.hfp,'userdata',ud);
    set(display.hfp,'visible','on');
    set(display.ht,'visible',within)
    set(display.hl,'visible',within)
    set(display.hbg,'visible',within)
    set(display.hr1,'visible',within)
    set(display.hr2,'visible',within)
    set(display.hb,'visible',within)
end

dim = options.dim;
if dim.n_phi > 0
    m = p_group.muPhi;
    v = diag(p_group.SigmaPhi);
    cla(display.ha(1,1))
    set(display.ha(1,1),'visible','on');
    [haf,hf,hp] = plotUncertainTimeSeries(m,v,1,display.ha(1,1));
    set(display.ha(1,1),'xlim',[0.5,options.dim.n_phi+0.5],'box','off');
    indrfx = setdiff(1:dim.n_phi,o_group.ind.phi_ffx);
    indrfx = intersect(indrfx,o_group.ind.phi_in);
    if ~isempty(indrfx)
        VBA_title(display.ha(1,1),'phi: population mean')
        cla(display.ha(1,2))
        set(display.ha(1,2),'visible','on');
        alphaHat = p_group.a_vPhi./p_group.b_vPhi;
        var_alpha = alphaHat./p_group.b_vPhi;
        logCI = log(alphaHat+sqrt(var_alpha)) - log(alphaHat);
        plotUncertainTimeSeries(log(alphaHat),logCI.^2,1,display.ha(1,2));
        set(display.ha(1,2),'xlim',[0.5,options.dim.n_phi+0.5],'box','off','xtick',[])
        VBA_title(display.ha(1,2),'phi: population log-precision')
    else
        VBA_title(display.ha(1,1),'phi: group mean (FFX)')
    end
end
if dim.n_theta > 0
    m = p_group.muTheta;
    v = diag(p_group.SigmaTheta);
    cla(display.ha(2,1))
    set(display.ha(2,1),'visible','on');
    [haf,hf,hp] = plotUncertainTimeSeries(m,v,1,display.ha(2,1));
    set(display.ha(2,1),'xlim',[0.5,options.dim.n_theta+0.5],'box','off');
    indrfx = setdiff(1:dim.n_theta,o_group.ind.theta_ffx);
    indrfx = intersect(indrfx,o_group.ind.theta_in);
    if ~isempty(indrfx)
        VBA_title(display.ha(2,1),'theta: population mean')
        cla(display.ha(2,2))
        set(display.ha(2,2),'visible','on');
        alphaHat = p_group.a_vTheta./p_group.b_vTheta;
        var_alpha = alphaHat./p_group.b_vTheta;
        logCI = log(alphaHat+sqrt(var_alpha)) - log(alphaHat);
        plotUncertainTimeSeries(log(alphaHat),logCI.^2,1,display.ha(2,2));
        set(display.ha(2,2),'xlim',[0.5,options.dim.n_theta+0.5],'box','off','xtick',[])
        VBA_title(display.ha(2,2),'theta: population log-precision')
    else
        VBA_title(display.ha(2,1),'theta: group mean (FFX)')
    end
end
if dim.n >0
    m = p_group.muX0;
    v = diag(p_group.SigmaX0);
    cla(display.ha(3,1))
    set(display.ha(3,1),'visible','on');
    [haf,hf,hp] = plotUncertainTimeSeries(m,v,1,display.ha(3,1));
    set(display.ha(3,1),'xlim',[0.5,options.dim.n+0.5],'box','off');
    indrfx = setdiff(1:dim.n,o_group.ind.x0_ffx);
    indrfx = intersect(indrfx,o_group.ind.x0_in);
    if ~isempty(indrfx)
        VBA_title(display.ha(3,1),'x0: population mean')
        cla(display.ha(3,2))
        set(display.ha(3,2),'visible','on');
        alphaHat = p_group.a_vX0./p_group.b_vX0;
        var_alpha = alphaHat./p_group.b_vX0;
        logCI = log(alphaHat+sqrt(var_alpha)) - log(alphaHat);
        plotUncertainTimeSeries(log(alphaHat),logCI.^2,1,display.ha(3,2));
        set(display.ha(3,2),'xlim',[0.5,options.dim.n+0.5],'box','off','xtick',[])
        VBA_title(display.ha(3,2),'x0: population log-precision')
    else
        VBA_title(display.ha(3,1),'x0: group mean (FFX)')
    end
end


% plot free energy update
if o_group.it >= 1
    cla(display.ha(4,1))
    set(display.ha(4,1),'visible','on');
    plot([0:o_group.it],o_group.F,'marker','.','parent',display.ha(4,1))
    VBA_title(display.ha(4,1),'VB optimization: F values')
    ylabel(display.ha(4,1),'free energy')
    xlabel(display.ha(4,1),'VB iterations')
    xl = {'prior(0)',['posterior(',num2str(o_group.it),')']};
    set(display.ha(4,1),'xtick',[0,o_group.it],'xticklabel',xl,'box','off');
end

VBA_getSubplots ();


function [] = VBAMFX_displayWithin(i1,i2)
try
    ud = get(get(i1,'parent'),'userdata');
    is = get(ud.display.hl,'value');
    if get(ud.display.hr1,'value')
        posterior = ud.p_sub{is};
        out = ud.o_sub{is};
        out.options.figName = ['VBA MFX: within-subject model inversion: subject #',num2str(is),' (empirical Bayes priors)'];
    else
        posterior = ud.o_group.initVBA.p_sub{is};
        out = ud.o_group.initVBA.o_sub{is};
        out.options.figName = ['VBA MFX: within-subject model inversion: subject #',num2str(is),' (initial priors)'];
    end
    [hfp,out] = VBA_ReDisplay(posterior,out,1);
catch
    disp(' --- Oops: something went wrong when displaying within-subject VBA model inversion! ---')
end


