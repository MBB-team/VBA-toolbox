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

% by default, show posterior statistics
if ~exist('priors','var')
    priors = false;
end

if ~options.DisplayWin
    return
else
    visible = 'on';
end

% check whether this is deterministic model
options0 = options;
if isequal(options.g_fname,@VBA_odeLim)
    options = options.inG.old.options;
end


%% set up figure and panel if necessary
% set up
if isfield(options,'display') 
    display = options.display ;
else
    display = struct;
end

% check if a figure has already been set up
if isfield(display,'hfp') % from init
    hfp = display.hfp;
elseif isfield(options,'hf') % form metaiteration
    hfp = options.hf ;
    try
        delete(get(hfp,'children'));
    end
else
    pos0 = get(0,'screenSize');
    pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];
    hfp = figure('position',pos,'color',[1 1 1],'menubar','none','tag','VBNLSS','Renderer','OpenGL');
end
display.hfp = hfp;
set(display.hfp,'name',options.figName);

% check if a panel has already been set up (spm tab or central element)
hPanel = getPanel(display.hfp);
if ~isempty(hPanel)
    % if there is one, start from scratch
    delete(get(hPanel,'children'));
else
    % otherwise, create a central panel to gather the plots
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
    xl = 'data dimension';
else
    xlim = [1,options.dim.n_t];
    xl = 'time';
end
if isequal(xlim,[1,1])
    xlim = [1-eps,1+eps];
end


% if multisource, show a select menu to change the displayed source
if numel(options.sources)>1
    % name the sources by their number
    for i=1:numel(options.sources)
        snames{i} = ['#',num2str(i)];
    end
    % show the select menu
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
        'visible'   ,visible                , ...
        'Value',currentSource               );
    % show a label text next to the menu
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
    %store handles
    display.source_selector = handles;
   
else
    display.source_selector = [];
end
    
%% prepare figure by initializing all axes
if priors
    data_label = 'prior' ;
    data_conditioner = 'm';
else
    data_label = 'posterior' ;
    data_conditioner = 'y,m';
end

% 1) predictive density and observations as function of time
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
h = stylePlot(subplot( ...
    'Position'  ,[.1 .75 .525 .165] , ...
    'parent'    ,hPanel             , ...
    'xlim'      ,xlim               , ...
    'ygrid'     ,'on'               , ...
    'xgrid'     ,'off'              ));
axis(h,'tight')

VBA_title(h, ...
    sprintf('%s predictive density: p(g(x)|%s)',data_label,data_conditioner)) ;
styleLabel(xlabel(h, ...
    sprintf('%s',xl)));
styleLabel(ylabel(h, ...
    sprintf('<g(x)|%s> & y',data_conditioner)));

display.ha(1) = h;

% 2) model fit, shown as prediciton against observation
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
h = stylePlot(subplot( ...
    'Position'  ,[.7 .75 .225 .165] , ...
    'parent'    ,hPanel             , ...
    'ygrid'     ,'on'               , ...
    'xgrid'     ,'on'               ));
axis(h,'tight')

VBA_title(h, ...
    sprintf('Model fit: <g(x)|%s> versus y',data_conditioner));
styleLabel(xlabel(h, ...
    sprintf('<g(x)|%s>',data_conditioner)));
styleLabel(ylabel(h, ...
    sprintf('y')));

display.ha(2) = h;

% 3) hidden states as function of time
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
h = stylePlot(subplot( ...
    'Position'  ,[.1 .51 .525 .165] , ...
    'parent'    ,hPanel             , ...
    'xlim'      ,xlim               , ...
    'ygrid'     ,'on'               , ...
    'xgrid'     ,'off'              ));
axis(h,'tight')
        
if options.dim.n == 0
    placeHolder(h,'no hidden states timeseries')
else
    VBA_title(h, ...
        sprintf('hidden states: p(x|%s)',data_conditioner));
    styleLabel(xlabel(h, ...
        sprintf('time')));
    styleLabel(ylabel(h, ...
        sprintf('<x | %s>',data_conditioner)));
end
if (options.dim.n_t == 1)
    styleLabel(xlabel(h, ...
        sprintf('state dimension')));
end
display.ha(3) = h;

% 4) Initial state
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
h = stylePlot(subplot( ...
    'Position'  ,[.7 .51 .225 .165]      , ...
    'parent'    ,hPanel                  , ...
    'xlim'      ,[0.2,options.dim.n+0.8] , ...
    'xtick'     ,[]                      , ...
    'ygrid'     ,'on'                    , ...
    'xgrid'     ,'off'                   ));
      
if options.dim.n == 0
    placeHolder(h,'no initial hidden states')
else
    VBA_title(h, ...
        sprintf('initial conditions: p(x_0|%s)',data_conditioner));

    if options.updateX0
        x_label_str = 'x_0 dimensions' ;
        if ~priors
            y_label_str = '<x_0 | y,m> - <x_0 | m>';
        else
            y_label_str = '<x_0 | m>' ;
        end
    else
        x_label_str = 'x_0 dimensions [fixed pdf]' ;
        y_label_str = 'x_0' ;
    end
    styleLabel(xlabel(h,x_label_str));
    styleLabel(ylabel(h,y_label_str));
    
end
display.ha(4) = h;

% 5) Observation parameters
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
h = stylePlot(subplot( ...
    'Position'  ,[.1 .27 .525 .165]          , ...
    'parent'    ,hPanel                      , ...
    'xlim'      ,[0.2,options.dim.n_phi+0.8] , ...
    'xtick'     ,[]                          , ...
    'ygrid'     ,'on'                        , ...
    'xgrid'     ,'off'                       ));

if options.dim.n_phi == 0
    placeHolder(h,'no observation parameters')
else
    VBA_title(h, ...
        sprintf('observation parameters: p(\\phi|%s)',data_conditioner));

    if ~options.OnLine
        x_label_str = '\phi dimensions' ;
    else
        x_label_str = 'time' ;
    end
    if ~priors
        y_label_str = '<\phi | y,m> - <\phi | m>' ;
    else
        y_label_str = '<\phi | m>' ;
    end
    styleLabel(xlabel(h,x_label_str));
    styleLabel(ylabel(h,y_label_str));
    
    plotUncertainTimeSeries(zeros(options.dim.n_phi,1),zeros(options.dim.n_phi,1),1,h);

end
display.ha(5) = h;

% 6) measurement noise precision hyperparameter
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Ngs=sum([options.sources(:).type]==0);
    
h = stylePlot(subplot( ...
    'Position'  ,[.7 .27 .225 .165] , ...
    'parent'    ,hPanel             , ...
    'xlim'      ,[0.2,Ngs+0.8]      , ...
    'xtick'     ,1:Ngs              , ...
    'ygrid'     ,'on'               , ...
    'xgrid'     ,'off'              ));

if Ngs == 0
    placeHolder(h,'no observation hyperparameters')
else
    VBA_title(h, ...
        sprintf('observation precision: p(\\sigma | %s)',data_conditioner));

        if ~options.OnLine && options.updateHP
            x_label_str = 'gaussian source' ;
        elseif ~options.OnLine && ~options.updateHP
            x_label_str = '[fixed pdf]' ;
        else
            x_label_str = 'time' ;
        end
        y_label_str = '<log(\sigma)>' ;
        styleLabel(xlabel(h,x_label_str));
        styleLabel(ylabel(h,y_label_str));
end
display.ha(6) = h;

% 7) Evolution parameters
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
h = stylePlot(subplot( ...
    'Position'  ,[.1 .03 .525 .165]             , ...
    'parent'    ,hPanel                         , ...
    'xlim'      ,[0.2,options.dim.n_theta+0.8]  , ...
    'xtick'     ,[]                             , ...
    'ygrid'     ,'on'                           , ...
    'xgrid'     ,'off'                          ));

if options.dim.n_theta == 0
    placeHolder(h,'no evolution parameters')
else
    VBA_title(h, ...
        sprintf('evolution parameters: p(\\theta | %s)',data_conditioner));

    if ~options.OnLine
        x_label_str = '\theta dimensions' ;
    else
        x_label_str = 'time' ;
    end
    if ~priors
        y_label_str = '<\theta | y,m> - <\theta | m>' ;
    else
        y_label_str = '<\theta | m>' ;
    end
    styleLabel(xlabel(h,x_label_str));
    styleLabel(ylabel(h,y_label_str));
    
end  
display.ha(7) = h;

% 8) State noise precision hyperparameter
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
h = stylePlot(subplot( ...
    'Position'  ,[.7 .03 .225 .165] , ...
    'parent'    ,hPanel             , ...
    'xlim'      ,[0.2 1.8]          , ...
    'xtick'     ,[]                 , ...
    'ygrid'     ,'on'               , ...
    'xgrid'     ,'off'              ));

if isequal(options0.g_fname,@VBA_odeLim) || options.dim.n == 0 || isinf(options0.priors.a_alpha/options0.priors.b_alpha) % not for non stochastic systems
    placeHolder(h,'no evolution hyperparameters')
else
    VBA_title(h, ...
        sprintf('evolution precision: p(\\alpha|%s)',data_conditioner));

    if ~options.OnLine && options.updateHP
        x_label_str = '' ;
    elseif ~options.OnLine && ~options.updateHP
        x_label_str = '[fixed pdf]' ;
    else
        x_label_str = 'time' ;
    end
    y_label_str = '<log(\alpha)>' ;
    styleLabel(xlabel(h,x_label_str));
    styleLabel(ylabel(h,y_label_str));
    
end 
display.ha(8) = h;

%% Create text boxes for user feedback
display.ho     = uicontrol('parent',display.hfp,'style','text','tag','VBLaplace','units','normalized','position',[0.20 0.010 0.60 0.02],'backgroundcolor',[1,1,1]);
display.hm(1)  = uicontrol('parent',display.hfp,'style','text','tag','VBLaplace','units','normalized','position',[0.28 0.035 0.40 0.02],'backgroundcolor',[1,1,1]);
display.hm(2)  = uicontrol('parent',display.hfp,'style','text','tag','VBLaplace','units','normalized','position',[0.68 0.035 0.10 0.02],'backgroundcolor',[1,1,1]);
display.htt(1) = uicontrol('parent',display.hfp,'style','text','tag','VBLaplace','units','normalized','position',[0.70 0.965 0.30 0.02],'backgroundcolor',[1,1,1]);

% Create 'pause' uicontrol button
if ~isfield(options,'noPause') || ~options.noPause
    vis = 'on';
else
    vis = 'off';
end
display.hpause = uicontrol('parent',display.hfp,'style','toggle','tag','VBLaplace','units','normalized','position',[0.4,0.96,0.2,0.02],'backgroundcolor',.8*[1,1,1],'string','pause and diagnose?','tag','pause_vb','visible',vis);

%% actually display
drawnow
try
    VBA_getSubplots ();
end

%% save handles and options
display.OnLine  = options.OnLine;
options         = options0;
options.display = display;


%% ========================================================================
% subfunctions

% callback for the source selector menu
function changeSource(hObject,evt,si)

    % get data stored in the figure
    hPanel = get(hObject,'parent');
    ud = get(hPanel,'userdata');
    
    % update the source index and store back
    ud.currentSource = get(hObject,'Value');
    set(hPanel,'userdata',ud);
    
    % try to update display
    if isfield(ud,'update_plot')
        ud.update_plot();
    end

% styling of axes labels
function h=stylePlot(h)
    set(h,'nextplot'  ,'add'         );
    set(h,'tag'       ,'VBLaplace'   );
    set(h,'box'       ,'off'         );
     
function styleLabel(h)
    set(h,'FontName'    ,'Arial'    );
    set(h,'FontUnits'   ,'points'   );
    set(h,'FontSize'    ,10         );
    set(h,'FontWeight'  ,'normal'   );
    set(h,'Color'       ,[0 0 0]    );
 
 % display low key text when usual plot is not required
 function placeHolder(h,label)
     xx = get(h,'XLim');
     yy = get(h,'YLim');
     t=text(mean(xx),mean(yy),label,'parent',h);
     set(t, ...
         'HorizontalAlignment','center'     , ...
         'FontSize'           ,10           , ...
         'Color'              ,[.6 .6 .6]   );
    set(h,'Visible','off');
     