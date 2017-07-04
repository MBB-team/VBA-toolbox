function [options] = VBA_initDisplay(options,priors)
% initialize VBA graphical output display
% function [options] = VBA_initDisplay(options,priors)
% This function is typically called when initializing VBA's graphical
% output display (e.g., at the begining of VBA inversion if
% options.DisplayWin=1). In this case, a subset of VBA posterior summary
% statistics are displayed. But it can also be called when VBA inversion is
% complete, for both posterior and prior summary statisitcs.
% IN:
%   - options: VBA's options structure [see VBA_NLStateSpaceModel.m]
%   - priors: flag for displaying prior summary statistics
% OUT:
%   - options: VBA_initDisplay inserts a 'display' field in VBA's options
%   structure. This field contains the handles of the relevant graphical
%   objects.

try
    priors;
catch
    priors = 0;
end

if ~options.DisplayWin
    return
else
    visible = 'on';
end

% First check whether this is standard DCM
options0 = options;
if isequal(options.g_fname,@VBA_odeLim)
    options = options.inG.old.options;
end


%% set up figure and panel if necessary
% set up
try 
    display = options.display ;
catch
    display = struct;
end

% check if a figure has already been set up
if isfield(display,'hfp') % from init
    hfp = display.hfp;
elseif isfield(options,'hf') % form metaiteration
    hfp = options.hf ;
    clo(hfp)
else
    pos0 = get(0,'screenSize');
    pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];
    hfp = figure('position',pos,'color',[1 1 1],'menubar','none','tag','VBNLSS','Renderer','OpenGL');
end
display.hfp = hfp;
set(display.hfp,'name',options.figName);

% check if a panel has already been set up (spm tab or central element)
hPanel = getPanel(display.hfp);
if isempty(hPanel)
    hPanel = uipanel('parent',display.hfp,'Tag','VBLaplace','BorderType','none','BackgroundColor',[1 1 1]);
    set(hPanel,'units','normalized');
    set(hPanel,'Position',[.02 .08 .96 .87]);
end

%check if source was selected
ud = get(hPanel,'userdata');
if ~isfield(ud,'currentSource')
    ud.currentSource = 1;
    set(hPanel,'userdata',ud);
end
currentSource = ud.currentSource;

%% Create axes

% Create axes for measured and predicted data
if options.dim.n_t == 1
    xlim = [1,options.dim.p];
    xl = 'data dimensions';
else
    xlim = [1,options.dim.n_t];
    xl = 'time';
end
if isequal(xlim,[1,1])
    xlim = [1-eps,1+eps];
end


% if multisource, show selector
if numel(options.sources)>1
    for i=1:numel(options.sources) %dim s
        snames{i} = ['#',num2str(i)];
    end
    handles(1) = uicontrol( ...
        'style'     ,'popupmenu'            , ...
        'parent'    ,hPanel                 , ...
        'tag'       ,'VBLaplace'            , ...
        'units'     ,'normalized'           , ...
        'position'  ,[.20 0.96 0.13 0.02]   , ...
        'FontSize'  ,12                     , ...
        'FontWeight','bold'                 , ...
        'string'    ,snames                 , ...
        'callback'  ,@changeSource          , ...
        'visible'   ,visible                );
    handles(2) = uicontrol( ...
        'style'     ,'text'                 , ...
        'parent'    ,hPanel                 , ...
        'tag'       ,'VBLaplace'            , ...
        'BackgroundColor' ,[1 1 1]          , ...
        'units'     ,'normalized'           , ...
        'position'  ,[0.05 0.96 0.15 0.02]  , ...
        'FontSize'  ,12                     , ...
        'FontWeight','bold'                 , ...
        'HorizontalAlignment','left'        , ...
        'string'    ,'observation source:'              , ...
        'visible'   ,visible                );
    
    set(handles(1),'Value',currentSource);

    display.source_selector = handles;
    
else
    display.source_selector = [];
end
    
display.ha(1) = subplot('Position',[.1 .75 .525 .165],'parent',hPanel,'xlim',xlim,'nextplot','add','tag','VBLaplace','box','off');
if ~priors
    title(display.ha(1),'posterior predictive density: p(g(x)|y,m)','fontsize',12)
else
    title(display.ha(1),'prior predictive density: p(g(x)|m)','fontsize',12)
end
xlabel(display.ha(1),xl,'fontsize',10)
if ~priors
    ylabel(display.ha(1),'<g(x)|y,m> & y','fontsize',10)
else
    ylabel(display.ha(1),'<g(x)|m> & y','fontsize',10)
end

display.ha(2) = subplot('Position',[.7 .75 .225 .165],'parent',hPanel,'nextplot','add','tag','VBLaplace','box','off');

if ~priors
    title(display.ha(2),'Model fit: <g(x)|y,m> versus y','fontsize',12)
    xlabel(display.ha(2),'<g(x)|y,m>','fontsize',10)
else
    title(display.ha(2),'Model fit: <g(x)|m> versus y','fontsize',12)
    xlabel(display.ha(2),'<g(x)|m>','fontsize',10)
end
ylabel(display.ha(2),'y','fontsize',10)

% Create axes for hidden states and initial conditions

if options.dim.n > 0
    display.ha(3) = subplot('Position',[.1 .51 .375 .165],'parent',hPanel,'xlim',xlim,'nextplot','add','tag','VBLaplace','box','off');

    if ~priors
        title(display.ha(3),'hidden states: p(x|y,m)','fontsize',12)
    else
        title(display.ha(3),'hidden states: p(x|m)','fontsize',12)
    end
    xlabel(display.ha(3),'time','fontsize',10)
    if ~priors
        ylabel(display.ha(3),'<x|y,m>','fontsize',10)
    else
        ylabel(display.ha(3),'<x|m>','fontsize',10)
    end
    display.ha(4) = subplot('Position',[.55 .51 .375 .165],'parent',hPanel,'nextplot','add','xlim',[0.2,options.dim.n+0.8],'xtick',[],'tag','VBLaplace','box','off');
    if ~priors
        title( display.ha(4),'initial conditions: p(x_0|y,m)','fontsize',12)
    else
        title( display.ha(4),'initial conditions: p(x_0|m)','fontsize',12)
    end
    if options.updateX0
        xlabel(display.ha(4),'x_0 dimensions','fontsize',10)
        if ~priors
            ylabel(display.ha(4),'<x_0|y,m> - <x_0|m>','fontsize',10)
        else
            ylabel(display.ha(4),'<x_0|m>','fontsize',10)
        end
    else
        xlabel(display.ha(4),'x_0 dimensions [fixed pdf]','fontsize',10)
        ylabel(display.ha(4),'x_0','fontsize',10)
    end
end

% Create axes for observation parameters
if options.dim.n_phi > 0
    display.ha(5) = subplot('Position',[.1 .27 .375 .165],'parent',hPanel,'nextplot','add','xlim',[0.2,options.dim.n_phi+0.8],'xtick',[],'tag','VBLaplace','box','off');

    if ~priors
        title(display.ha(5),'observation parameters: p(phi|y,m)','fontsize',12)
    else
        title(display.ha(5),'observation parameters: p(phi|m)','fontsize',12)
    end
    if ~options.OnLine
        xlabel(display.ha(5),'phi dimensions','fontsize',10)
    else
        xlabel(display.ha(5),'time','fontsize',10)
    end
    if ~priors
        ylabel(display.ha(5),'<phi|y,m> - <phi|m>','fontsize',10)
    else
        ylabel(display.ha(5),'<phi|m>','fontsize',10)
    end
end

% Create axes for measurement noise precision hyperparameter
Ngs=sum([options.sources(:).type]==0);
if Ngs>0
    display.ha(6) = subplot('Position',[.55 .27 .375 .165],'parent',hPanel,'nextplot','add','xlim',[0.2,Ngs+0.8],'xtick',1:Ngs,'tag','VBLaplace','box','off');
    if ~priors
        title(display.ha(6),'measurement noise precision: p(sigma|y,m)','fontsize',12)
    else
        title(display.ha(6),'measurement noise precision: p(sigma|m)','fontsize',12)
    end
    if ~options.OnLine && options.updateHP
        xlabel(display.ha(6),'gaussian source','fontsize',10)
    elseif ~options.OnLine && ~options.updateHP
        xlabel(display.ha(6),'[fixed pdf]','fontsize',10)
    else
        xlabel(display.ha(6),'time','fontsize',10)
    end
    ylabel(display.ha(6),'<log(sigma)>','fontsize',10)
end

% Create axes for evolution parameters
if options.dim.n_theta > 0
    display.ha(7) = subplot('Position',[.1 .03 .375 .165],'parent',hPanel,'nextplot','add','xlim',[0.2,options.dim.n_theta+0.8],'xtick',[],'tag','VBLaplace','box','off');
    
    if ~priors
        title(display.ha(7),'evolution parameters: p(theta|y,m)','fontsize',12)
    else
        title(display.ha(7),'evolution parameters: p(theta|m)','fontsize',12)
    end
    if ~options.OnLine
        xlabel(display.ha(7),'theta dimensions','fontsize',10)
    else
        xlabel(display.ha(7),'time','fontsize',10)
    end
    if ~priors
        ylabel(display.ha(7),'<theta|y,m> - <theta|m>','fontsize',10)
    else
        ylabel(display.ha(7),'<theta|m>','fontsize',10)
    end
end

% Create axes for state noise precision hyperparameter
if ~isequal(options0.g_fname,@VBA_odeLim) && options.dim.n > 0 % not for non stochastic systems
    display.ha(8) = subplot('Position',[.55 .03 .375 .165],'parent',hPanel,'xlim',[0.2,1.8],'xtick',[],'nextplot','add','tag','VBLaplace','box','off');
    if ~priors
        title(display.ha(8),'system''s noise precision: p(alpha|y,m)','fontsize',12)
    else
        title(display.ha(8),'system''s noise precision: p(alpha|m)','fontsize',12)
    end
    if ~options.OnLine && options.updateHP
        xlabel(display.ha(8),'','fontsize',10)
    elseif ~options.OnLine && ~options.updateHP
        xlabel(display.ha(8),'[fixed pdf]','fontsize',10)
    else
        xlabel(display.ha(8),'time','fontsize',10)
    end
    ylabel(display.ha(8),'<log(alpha)>','fontsize',10)
end

% Create text boxes for user feedback
display.ho = uicontrol('parent',display.hfp,'style','text','tag','VBLaplace','units','normalized','position',[0.2,0.01,0.6,0.02],'backgroundcolor',[1,1,1]);
display.hm(1) = uicontrol('parent',display.hfp,'style','text','tag','VBLaplace','units','normalized','position',[0.28,0.035,0.4,0.02],'backgroundcolor',[1,1,1]);
display.hm(2) = uicontrol('parent',display.hfp,'style','text','tag','VBLaplace','units','normalized','position',[0.68,0.035,0.1,0.02],'backgroundcolor',[1,1,1]);
display.htt(1) = uicontrol('parent',display.hfp,'style','text','tag','VBLaplace','units','normalized','position',[0.75 0.97 0.25 0.02],'backgroundcolor',[1,1,1]);

% Create 'pause' uicontrol button
try
    if ~options.noPause
        vis = 'on';
    else
        vis = 'off';
    end
catch
    vis = 'on';
end
display.hpause = uicontrol('parent',display.hfp,'style','toggle','tag','VBLaplace','units','normalized','position',[0.4,0.96,0.2,0.02],'backgroundcolor',.8*[1,1,1],'string','pause and diagnose?','tag','pause_vb','visible',vis);


%set(hPanel,'userdata',currentSource);

drawnow
display.OnLine = options.OnLine;
options = options0;
options.display = display;



try
    getSubplots
end


function changeSource(hObject,evt,si)

    hPanel = get(hObject,'parent');
    ud = get(hPanel,'userdata');
    
    ud.currentSource = get(hObject,'Value');
    set(hPanel,'userdata',ud);
    
    if isfield(ud,'update_plot')
        ud.update_plot();
    end
