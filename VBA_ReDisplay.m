function [hfp,out] = VBA_ReDisplay(posterior,out,newFig,fromPause)
% re-creates the graphical output of the VBA inversion + diagnostics
% function [hfp] = VBA_ReDisplay(posterior,out,newFig)
% VBA_ReDisplay first looks for a figure with a tag 'VBNLSS', i.e. a figure
% that was already opened to review a model inversion, and clears it if it
% finds it (except if newFig=1).
% IN:
%   - posterior/out: standard output of VBA_NLStateSpaceModel.m
%   - newFig: a flag for creating a new VBA display figure.
% OUT:
%   - hfp: handle of the display figure
%   - out: standard output of VBA_NLStateSpaceModel.m, augmented with
%   diagnostics, if those were not included in the out structure before the
%   call to VBA_ReDisplay.m.

% set default
if ~exist('newFig','var')
    newFig=0;
end
if ~exist('fromPause','var')
    fromPause=0;
end

% create shortcuts
options = out.options;

% check if a figure already exist or create a new one if needed or requested
hfp = findobj('tag','VBNLSS');
if isempty(hfp) || newFig
    pos0 = get(0,'screenSize');
    pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];
    hfp = figure( ...
        'position'  , pos               , ...
        'color'     , [1 1 1]           , ...
        'name'      , options.figName   , ...
        'menubar'   , 'none'            , ...
        'tag'       , 'VBNLSS'          , ...
        'Renderer'  , 'OpenGL'          );
else
    hfp = hfp(1);
    clf(hfp)
    set(hfp,'name',options.figName);
end


% ensure display in switched on
options.DisplayWin = 1;



% compute missing statistics if needed
if ~isfield(out,'diagnostics')
    out.diagnostics = VBA_getDiagnostics(posterior,out);
end

% store data in figure for later 
ud.posterior = posterior;
ud.out = out;

set(hfp,'userdata',ud);

% setup tabs

if ~isempty(out.diagnostics.kernels)
    labels = {'summary','VB inversion','diagnostics','kernels','conv','priors'};
    callbacks = {@mySummary,@myVB,@myDiagnostics,@myKernels,@myConv,@myPriors};
else
    labels = {'summary','VB inversion','diagnostics','conv','priors'};
    callbacks = {@mySummary,@myVB,@myDiagnostics,@myConv,@myPriors};
end
if out.dim.n > 0 && ~isinf(out.options.priors.a_alpha) && ~isequal(out.options.priors.b_alpha,0) && ~out.options.OnLine
    labels{end+1} = 'deterministic';
    callbacks{end+1} = @myDeterministic;
end

if fromPause
    active = 2;
else
    active = 1;
end

[handles] = VBA_spm_uitab(hfp,labels,callbacks,'diagnostics_tabs',active);
set(handles.htab   , 'backgroundcolor' , [1 1 1]     );
set(handles.hh     , 'backgroundcolor' , [1 1 1]     );
set(handles.hp     , 'backgroundcolor' , [1 1 1]     );
set(handles.hp     , 'HighlightColor'  , 0.8*[1 1 1] );
set(handles.htab(1), 'tooltipstring'   , 'summary description of the VB inversion' );
set(handles.htab(2), 'tooltipstring'   , 'results of the VB inversion (posterior pdfs)');
set(handles.htab(3), 'tooltipstring'   , 'VB inversion diagnostics (residuals and parameters covariance matrices)');
if ~isempty(out.diagnostics.kernels)
    ind = 4;
    set(handles.htab(4),'tooltipstring','system''s 1st-order Volterra kernels')
else
    ind = 3;
end
set(handles.htab(ind+1),'tooltipstring','history of free energy values along VB optimization')
set(handles.htab(ind+2),'tooltipstring','priors and associated predictive densities (under the Laplace assumption)')
if out.dim.n > 0 && ~isinf(out.options.priors.a_alpha) && ~isequal(out.options.priors.b_alpha,0) && ~out.options.OnLine
    set(handles.htab(ind+3),'tooltipstring','results of the VB inversion of the deterministic system')
end

if fromPause
    feval(@myVB,hfp)
else
    feval(@mySummary,hfp)
end

% =========================================================================
% tabs callbacks
% =========================================================================

function mySummary(hfig) 

    try 
        hfig;
    catch
        hfig = get(gco,'parent');
    end
    cleanPanel(hfig);

    ud = get(hfig,'userdata');
    out = ud.out;

    str = VBA_summary(out,1);
    for i=1:length(str)
        str{i} = sprintf(str{i});
    end
    str{7} = sprintf(['Estimation efficiency (minus posterior entropies):','\n ']);
    if ~isnan(out.diagnostics.efficiency.X)
        str{7} = sprintf([str{7},'    - hidden states: ',num2str(out.diagnostics.efficiency.X,'%4.3e'),'\n ']);
    end
    if ~isnan(out.diagnostics.efficiency.X0)
        str{7} = sprintf([str{7},'    - initial conditions: ',num2str(out.diagnostics.efficiency.X0,'%4.3e'),'\n ']);
    end
    if ~isnan(out.diagnostics.efficiency.Theta)
        str{7} = sprintf([str{7},'    - evolution parameters: ',num2str(out.diagnostics.efficiency.Theta,'%4.3e'),'\n ']);
    end
    if ~isnan(out.diagnostics.efficiency.Phi)
        str{7} = sprintf([str{7},'    - observation parameters: ',num2str(out.diagnostics.efficiency.Phi,'%4.3e'),'\n ']);
    end
    if ~isnan(out.diagnostics.efficiency.alpha)
        str{7} = sprintf([str{7},'    - state noise precision hyperparameter: ',num2str(out.diagnostics.efficiency.alpha,'%4.3e'),'\n ']);
    end
    if ~any(isnan(out.diagnostics.efficiency.sigma))
        gsi = find([out.options.sources.type]==0);
        sig_str = catnum2str(out.diagnostics.efficiency.sigma,gsi,length(gsi)>1);
        str{7} = sprintf([str{7},'    - data noise precision hyperparameter: ',sig_str,'\n ']);
    end
    uicontrol( ...
        'parent'                , hfig                  , ...
        'style'                 , 'text'                , ... 
        'tag'                   , 'VBLaplace'           , ... 
        'units'                 , 'normalized'          , ... 
        'position'              , [0.1,0.05,0.8,0.85]   , ... 
        'backgroundcolor'       , [1,1,1]               , ... 
        'HorizontalAlignment'   , 'left'                , ... 
        'fontsize'              , 11                    , ... 
        'string'                , str                   );

function myDeterministic(hfig)

    try hfig; catch, hfig = get(gco,'parent'); end
    cleanPanel(hfig);
    ud = get(hfig,'userdata');

    % Second: re-display VB-Laplace inversion output

    y         = ud.out.y;
    posterior = ud.out.options.init.posterior;
    out       = ud.out.options.init.out;

    options            = out.options;
    options.noPause    = 1;
    options.DisplayWin = 1;
    dim                = out.dim;
    suffStat           = out.suffStat;
    posterior.a_alpha  = Inf;
    posterior.b_alpha  = 0;

    % Initialize display figure
    options.display.hfp = hfig;
    options.figName     = get(hfig,'name');
    [options] = VBA_initDisplay(options);
    delete(options.display.htt)
    delete(options.display.hpause)
    delete(options.display.hm)
    delete(options.display.ho)
  
    hfig = options.display.hfp;
    drawnow

    % Display 
    if options.dim.n > 0
        VBA_updateDisplay(posterior,suffStat,options,y,0)
    end
    
    try, VBA_getSubplots (); end

    uicontrol( ...
        'parent'          , hfig                        , ...
        'style'           , 'pushbutton'                , ...
        'tag'             , 'VBLaplace'                 , ...
        'units'           , 'normalized'                , ...
        'position'        , [0.40 0.93 0.20 0.02]       , ...
        'backgroundcolor' , .8*[1,1,1]                  , ...
        'string'          , 'diagnose deterministic?'	, ...
        'callback'        , @diagnoseDeterministic      );

function diagnoseDeterministic(ho,e)
    ud = get(get(ho,'parent'),'userdata');
    posterior = ud.out.options.init.posterior;
    out = ud.out.options.init.out;
    VBA_ReDisplay(posterior,out,1);

function myPriors(hfig)

    try hfig; catch, hfig = get(gco,'parent'); end
    cleanPanel(hfig);
    ud = get(hfig,'userdata');

    % Second: re-display VB-Laplace inversion output
    out         = ud.out              ;
    y           = out.y               ;
    posterior   = out.options.priors  ;
    options     = out.options         ;
    options.noPause    = 1            ;
    options.DisplayWin = 1            ;
    dim         = out.dim             ;
    suffStat    = out.suffStat        ;
    suffStat.gx = out.diagnostics.pgx ;

    % set dx = -prior.muX (for display purposes)
    suffStat.dx0    = -posterior.muX0       ;
    suffStat.dtheta = -posterior.muTheta    ;
    suffStat.dphi   = -posterior.muPhi      ;
    suffStat.vy     =  out.diagnostics.pvy  ;

    % Initialize display figure
    options.display.hfp = hfig;
    options.figName     = get(hfig,'name');
    [options] = VBA_initDisplay(options,1);
    delete(options.display.htt)
    delete(options.display.hpause)
    delete(options.display.hm)
    delete(options.display.ho)
 
    % Display 
    options.OnLine = 0;
    VBA_updateDisplay(posterior,suffStat,options,y,0)

    try, VBA_getSubplots (); end

function myConv(hfig)

    try hfig; catch, hfig = get(gco,'parent'); end
    cleanPanel(hfig);
    ud     = get(hfig,'userdata');
    hPanel = getPanel(hfig);

    out         = ud.out;
    diagnostics = out.diagnostics;

    if length(out.suffStat.F)>2
        nit = length(out.suffStat.F)-1;
        ha = axes('parent',hPanel,'units','normalized','tag','VBLaplace','position',[0.15,0.6,0.5,0.3],'nextplot','add','xlim',[0,nit],'xtick',[0,nit],'xticklabel',{'prior (0)',['posterior (',num2str(nit),')']},'xgrid','off','ygrid','on');
        plot(ha,[0:nit],out.suffStat.F) 
        plot(ha,[0:nit],out.suffStat.F,'.') ;
        [haf,hp1,hp2] = plotUncertainTimeSeries(diagnostics.LLH0*ones(1,2),3^2*ones(1,2),[0,nit],ha);
        set(hp1,'facecolor',[1 0 0])
        set(hp2,'color',[1 0 0])
        text(nit/2,diagnostics.LLH0-3/2,'log p(y|H0)','color',[1 0 0],'parent',ha);
        if ~out.options.OnLine
            VBA_title(ha,'VB optimization: F values')
            xlabel(ha,'inner (Gauss-Newton) iterations')
        else
            VBA_title(ha,'online VB: F values')
            xlabel(ha,'time samples')
        end
        ylabel(ha,'Free energy')
        box(ha,'off')
        if ~out.options.OnLine
            xtl = {'first iteration','last iteration'};
        else
            xtl = {'first time point','last time point'};
        end
        ha = axes('parent',hPanel,'units','normalized','tag','VBLaplace','position',[0.15,0.15,0.5,0.3],'nextplot','add','xlim',[0,nit-1],'xtick',[0,nit-1],'xticklabel',xtl,'xgrid','off','ygrid','on');
        plot(ha,[0:nit-1],diff(out.suffStat.F)) ;
        plot(ha,[0:nit-1],diff(out.suffStat.F),'.') ;
        if ~out.options.OnLine
            VBA_title(ha,'VB optimization: F increments')
            xlabel(ha,'inner (Gauss-Newton) iterations')
        else
            VBA_title(ha,'online VB: F increments')
            xlabel(ha,'time samples')
        end
        ylabel(ha,'Free energy differences')
        box(ha,'off')
        try
            VBA_getSubplots ();
        end
    end

    options = orderfields(out.options);
    options = rmfield(options,'tStart');
    options = rmfield(options,'checkGrads');
    options = rmfield(options,'verbose');
    options = rmfield(options,'delays');
    options = rmfield(options,'isYout');
    options = rmfield(options,'skipf');
    finames = fieldnames(options);
    nopt = length(finames);
    str = {'Optional fields:';' '};
    for i=1:nopt
        tmp = getfield(options,finames{i});
        if ~isempty(tmp) && isnumeric(tmp)
            str{end+1} = [finames{i},' = ',num2str(max(tmp))];
        end
    end
    uicontrol('parent',hfig,'style','text','tag','VBLaplace','units','normalized','position',[0.75,0.1,0.2,0.8],'backgroundcolor',[1,1,1],'HorizontalAlignment','left','fontsize',11,'string',str);



function myKernels(hfig)
    
    try hfig; catch, hfig = get(gco,'parent'); end
    cleanPanel(hfig);
    ud     = get(hfig,'userdata');

    du = size(ud.out.diagnostics.kernels.y.m,3);
    unames = cell(du,1);
    for i=1:du %dim u
        unames{i} = ['#',num2str(i)];
    end
    handles(1) = uicontrol('style','popupmenu','parent',hfig,'tag','VBLaplace','units','normalized','position',[0.85 0.9 0.10 0.02],'fontsize',12,'string',unames,'callback',@myKerneli);
    handles(2) = uicontrol('style','text','parent',hfig,'tag','VBLaplace','BackgroundColor',get(hfig,'color'),'units','normalized','position',[0.82 0.93 0.16 0.02],'fontsize',12,'string','display input...');
    feval(@myKerneli,handles(1),[])



function myKerneli(hObject,evt)
    hfig = get(hObject,'parent');
    ind = get(hObject,'Value');
    ud = get(hfig,'userdata');
    try
        if isequal(get(ud.handles.hkernels(2),'parent'),hfig)
            delete(ud.handles.hkernels) ;
        end
    end
    out = ud.out;
    kernels = out.diagnostics.kernels;

    hPanel = getPanel(hfig);
    if ~isempty(kernels.x)
        % input effects - hidden states
        handles.hkernels(1) = subplot(2,1,1,'parent',hPanel,'nextplot','add','ygrid','on','tag','VBLaplace');
        pos = get(handles.hkernels(1),'position');
        set(handles.hkernels(1),'position',[0.2 pos(2) 0.6 pos(4)]);
        [t1,t2,hp] = plotUncertainTimeSeries(kernels.x.m(:,:,ind),kernels.x.v(:,:,ind),[],handles.hkernels(1));
        set(hp,'marker','.')
        set(handles.hkernels(1),'XLim',[0.5 size(kernels.x.m,2)+0.5],'xtick',[1:size(kernels.x.m,2)],'xticklabel',[0:size(kernels.x.m,2)-1])
        VBA_title(handles.hkernels(1),['states'' Volterra kernels: input #',num2str(ind),' (R2=',num2str(mean(kernels.x.R2),'%4.2f'),')'])
        ylabel(handles.hkernels(1),'(lagged) input weight')
        xlabel(handles.hkernels(1),'time lag')
    end
    % input effects - observables
    handles.hkernels(2) = subplot(2,1,2,'parent',hPanel,'nextplot','add','ygrid','on','tag','VBLaplace');
    pos = get(handles.hkernels(2),'position');
    set(handles.hkernels(2),'position',[0.2 pos(2) 0.6 pos(4)])
    hold(handles.hkernels(2),'on')
    plot(handles.hkernels(2),kernels.g.m(1,:,ind)','marker','.','color',[0 0 0])
    plot(handles.hkernels(2),kernels.y.m(1,:,ind)','marker','.','linestyle',':','color',[0 0 0])
    legend(handles.hkernels(2),{['simulated observables',' (R2=',num2str(mean(kernels.g.R2),'%4.2f'),')'],['observed samples',' (R2=',num2str(mean(kernels.y.R2),'%4.2f'),')']})
    plot(handles.hkernels(2),kernels.y.m(:,:,ind)','marker','.','linestyle',':')
    [t1,t2,hp] = plotUncertainTimeSeries(kernels.g.m(:,:,ind),kernels.g.v(:,:,ind),[],handles.hkernels(2));
    set(hp,'marker','.') ;
    set(handles.hkernels(2),'XLim',[0.5 size(kernels.g.m,2)+0.5],'xtick',[1:size(kernels.g.m,2)],'xticklabel',[0:size(kernels.g.m,2)-1]) ;
    VBA_title(handles.hkernels(2),['observables'' Volterra kernels: input #',num2str(ind)])
    ylabel(handles.hkernels(2),'(lagged) input weight')
    xlabel(handles.hkernels(2),'time lag')

    ud.handles = handles;
    set(hfig,'userdata',ud);
    try, VBA_getSubplots (); end

function myDiagnostics(hfig)

    try hfig; catch, hfig = get(gco,'parent'); end
    cleanPanel(hfig);
    ud     = get(hfig,'userdata');

    out         = ud.out;
    y           = out.y;
    diagnostics = out.diagnostics;

    hPanel = getPanel(hfig);
    % display micro-time hidden-states
    if ~isempty(diagnostics.MT_x)
        display.ha(1) = subplot(4,2,1,'parent',hPanel,'nextplot','add','tag','VBLaplace','ygrid','on','box','off');
        VBA_title(display.ha(1),'micro-time resolution predicted data')
        xlabel(display.ha(1),'time','fontsize',8)
        ylabel(display.ha(1),'g(x) & y','fontsize',8)
        plot(display.ha(1),diagnostics.microTime,diagnostics.MT_gx')
        plot(display.ha(1),diagnostics.microTime(diagnostics.sampleInd),diagnostics.MT_gx(:,diagnostics.sampleInd)','.')
        plot(display.ha(1),diagnostics.microTime(diagnostics.sampleInd),y,':')
        axis(display.ha(1),'tight')
        display.ha(2) = subplot(4,2,2,'parent',hPanel,'nextplot','add','tag','VBLaplace','ygrid','on','box','off');
        VBA_title(display.ha(2),'micro-time resolution hidden states')
        xlabel(display.ha(2),'time','fontsize',8)
        ylabel(display.ha(2),'x','fontsize',8)
        plot(display.ha(2),diagnostics.microTime,diagnostics.MT_x')
        plot(display.ha(2),diagnostics.microTime(diagnostics.sampleInd),diagnostics.MT_x(:,diagnostics.sampleInd)','.')
        axis(display.ha(2),'tight')
    end

    if numel(out.options.sources)>1
        vis = 'on';
    else
        vis = 'off';
    end
    ds = numel(ud.out.diagnostics.dy);
    snames = cell(ds,1);
    for i=1:ds %dim s
        snames{i} = ['#',num2str(i)];
    end
    handles(1) = uicontrol('style','popupmenu','parent',hfig,'tag','VBLaplace','units','normalized','position',[0.55 0.5 0.10 0.02],'fontsize',12,'string',snames,'callback',@myDiagnosticsi,'visible',vis);
    handles(2) = uicontrol('style','text','parent',hfig,'tag','VBLaplace','BackgroundColor',get(hfig,'color'),'units','normalized','position',[0.52 0.53 0.16 0.02],'fontsize',12,'string','source:','visible',vis);
    feval(@myDiagnosticsi,handles(1),[])

    % display state noise
    if ~isempty(diagnostics.dx.dx)
        xlim = [diagnostics.dx.nx(1)-diagnostics.dx.d,diagnostics.dx.nx(end)+diagnostics.dx.d];
        display.ha(4) = subplot(4,2,6,'parent',hPanel,'nextplot','add','xlim',xlim,'ygrid','on','tag','VBLaplace','box','off');
        VBA_title(display.ha(4),'weighted state noise distribution')
        xlabel(display.ha(4),'eta(t) = x(t+1)-f(x(t))','fontsize',8)
        ylabel(display.ha(4),'p(eta|y)','fontsize',8)
        bar(diagnostics.dx.nx,diagnostics.dx.ny,'facecolor',[.8 .8 .8],'parent',display.ha(4))
        plot(display.ha(4),diagnostics.dx.grid,diagnostics.dx.pg,'r')
        plot(display.ha(4),diagnostics.dx.grid,diagnostics.dx.pg2,'g')
        legend(display.ha(4),{'empirical histogram','Gaussian approx','posterior approx'})
        display.ha(8) = subplot(4,2,4,'parent',hPanel,'nextplot','add','tag','VBLaplace','ygrid','on','box','off');
        try
            plotUncertainTimeSeries(out.suffStat.dx,out.suffStat.vdx,diagnostics.microTime(diagnostics.sampleInd),display.ha(8));
        catch
            plot(display.ha(8),diagnostics.microTime(diagnostics.sampleInd),out.suffStat.dx','marker','.')
        end
        axis(display.ha(8),'tight')
        VBA_title(display.ha(8),'state noise time series')
        xlabel(display.ha(8),'time','fontsize',8)
        ylabel(display.ha(8),'eta(t) = x(t+1)-f(x(t))','fontsize',8)
    end

    % display parameters posterior correlation matrix
    display.ha(6) = subplot(4,2,8,'parent',hPanel);
    imagesc(diagnostics.C,'parent',display.ha(6))
    VBA_title(display.ha(6),'parameters posterior correlation matrix')
    set(display.ha(6),'tag','VBLaplace','xtick',diagnostics.ltick,'ytick',diagnostics.ltick,'xticklabel',diagnostics.ticklabel,'yticklabel',diagnostics.ticklabel,'box','off','nextplot','add');
    for i=1:length(diagnostics.tick)
        plot(display.ha(6),[0.5 size(diagnostics.C,1)+0.5],[diagnostics.tick(i) diagnostics.tick(i)],'color',[1 1 1])
        plot(display.ha(6),[diagnostics.tick(i) diagnostics.tick(i)],[0.5 size(diagnostics.C,1)+0.5],'color',[1 1 1])
    end
    grid(display.ha(6),'off')
    axis(display.ha(6),'square')
    set(display.ha(6),'clim',[-34/32 1]);
    col = colormap('jet');
    col(1,:) = 0.5*ones(1,3);
    colormap(display.ha(6),col);
    try display.hc(2) = colorbar('peer',display.ha(6)); end

    try VBA_getSubplots (); end

    
function myDiagnosticsi(hObject,evt,si)
    hfig = get(hObject,'parent');
    ind = get(hObject,'Value');
    ud = get(hfig,'userdata');
    hPanel = getPanel(hfig);
    try
        if isequal(get(ud.handles.hdiagnostics(1),'parent'),hPanel)
            delete(ud.handles.hdiagnostics)
        end
    end
    out = ud.out;
    dy = out.diagnostics.dy(ind);
    sourceType = out.options.sources(ind).type;
    sourceYidx = out.options.sources(ind).out;
    % display data noise
    xlim = [dy.nx(1)-dy.d,dy.nx(end)+dy.d];
    handles.hdiagnostics(1) = subplot(4,2,5,'parent',hPanel,'nextplot','add','xlim',xlim,'ygrid','on','tag','VBLaplace','box','off');
    VBA_title(handles.hdiagnostics(1),'weighted residuals distribution')
    xlabel(handles.hdiagnostics(1),'e(t) = y(t)-g(x(t))','fontsize',8)
    ylabel(handles.hdiagnostics(1),'p(e|y)','fontsize',8)
    bar(dy.nx,dy.ny,'facecolor',[.8 .8 .8],'parent',handles.hdiagnostics(1))
    plot(handles.hdiagnostics(1),dy.grid,dy.pg,'r')
    if sourceType==0
        plot(handles.hdiagnostics(1),dy.grid,dy.pg2,'g')
    end
    if sourceType==0
        legend(handles.hdiagnostics(1),{'empirical histogram','Gaussian approx','posterior approx'})
    else
        legend(handles.hdiagnostics(1),{'empirical histogram','Gaussian approx'})
    end
    if out.options.dim.n > 0
        gri = out.diagnostics.microTime(out.diagnostics.sampleInd);
        ti = 'time';
    else
        if out.options.dim.n_t>1
            gri = 1:out.options.dim.n_t;
            ti = 'time';
        else
            gri = 1:out.options.dim.p;
            ti = 'data dimensions';
        end
    end
    handles.hdiagnostics(2) = subplot(4,2,3,'parent',hPanel,'nextplot','add','tag','VBLaplace','ygrid','on','box','off');
    plot(handles.hdiagnostics(2),gri,out.suffStat.dy(sourceYidx,:)','marker','.')
    axis(handles.hdiagnostics(2),'tight')
    VBA_title(handles.hdiagnostics(2),'residuals time series')
    xlabel(handles.hdiagnostics(2),ti,'fontsize',8)
    ylabel(handles.hdiagnostics(2),'e(t) = y(t)-g(x(t))','fontsize',8)
    % display autocorrelation of residuals
    if ~ VBA_isWeird (dy.R) && out.dim.n_t > 1
        handles.hdiagnostics(3) = subplot(4,2,7,'parent',hPanel);
        plot(handles.hdiagnostics(3),[-out.options.dim.n_t:out.options.dim.n_t-1],fftshift(dy.R)')
        axis(handles.hdiagnostics(3),'tight')
        VBA_title(handles.hdiagnostics(3),'residuals empirical autocorrelation')
        xlabel(handles.hdiagnostics(3),'lag tau','fontsize',8)
        ylabel(handles.hdiagnostics(3),'Corr[e(t),e(t+tau)]','fontsize',8)
        set(handles.hdiagnostics(3),'tag','VBLaplace','ygrid','on','box','off');
    end
    ud.handles = handles;
    set(hfig,'userdata',ud);
    try, VBA_getSubplots (); end


function myVB(hfig)

    try hfig; catch, hfig = get(gco,'parent'); end
    cleanPanel(hfig);
    ud     = get(hfig,'userdata');

    out = ud.out;
    y   = out.y;
    posterior = ud.posterior;
    options = out.options;
    options.noPause = 1;
    options.DisplayWin =1;
    suffStat = out.suffStat;

    dim = out.dim;
    % Initialize display figure
    options.display.hfp = hfig;
    options.figName = get(hfig,'name');
    [options] = VBA_initDisplay(options);
    delete(options.display.htt)
    delete(options.display.hpause)
    delete(options.display.hm)
    delete(options.display.ho)
  
    hfig = options.display.hfp;
    drawnow
    
    % Display 
    VBA_updateDisplay(posterior,suffStat,options,y,0)

    try VBA_getSubplots (); end


%% helpers
function cleanPanel(hfig)
    hc = intersect(findobj('tag','VBLaplace'),get(hfig,'children'));
    if ~isempty(hc)
        delete(hc)
    end
    hPanel = getPanel(hfig);
    hc = intersect(findobj('tag','VBLaplace'),get(hPanel,'children'));
    if ~isempty(hc)
        delete(hc)
    end

function str = catnum2str(x,ind,many)
    str = [];
    for i=1:length(ind)
        si=ind(i);
        str = [str,', ',num2str(x(si),'%4.3e')];
        if many
            str = [str,' (source #',num2str(si),')'];
        end
    end
    str(1:2) = [];

