function VBA_updateDisplay(posterior,suffStat,options,y,it,flag)
% updates display of sufficient statistics
% function VBA_updateDisplay(F,posterior,suffStat,options,y,it,display,flag)
% This function deals with the screen display of iterative sufficient
% statistics updates of the VBA inversion algorithm

if ~options.DisplayWin
    return
end
display = options.display;

F = real(suffStat.F);

% check whether 'pause' button is toggled on
VBA_pause(options)

% replace source selector callback if needed
ud = get(getPanel(display.hfp),'userdata') ;

    function update_plot ()
        cla(display.ha(1));
        cla(display.ha(2));
        VBA_updateDisplay(posterior,suffStat,options,y,0,'Y');
    end

ud.update_plot = @update_plot; 
set(getPanel(display.hfp),'userdata',ud) ;

ud = VBA_check_struct(ud,'currentSource', 1 );

currentSource = ud.currentSource;

% 0- Check whether this is a deterministic dynamical system inversion
if isequal(options.g_fname,@VBA_odeLim)
    
    % Rebuild posterior from dummy 'ODE' posterior
    options0 = options;
    [posterior,options,dim,suffStat] = VBA_odeLim2NLSS(posterior,options,options.dim,suffStat,[],0);
    options.display = options0.display;
    
    % Then call VBA_updateDisplay again
    if ~isempty(it)
        VBA_updateDisplay(posterior,suffStat,options,y,it,'precisions')
    end
    VBA_updateDisplay(posterior,suffStat,options,y,it,'phi')
    VBA_updateDisplay(posterior,suffStat,options,y,it,'X')
    VBA_updateDisplay(posterior,suffStat,options,y,it,'theta')
    
    return
    
end


% Get sufficient statistics to be displayed
dTime = 1:size(y,2);
try
    gx = suffStat.gx(:,dTime);
    vy = suffStat.vy(:,dTime);
catch
    gx = nan(options.dim.p,numel(dTime));
    vy = nan(options.dim.p,numel(dTime));
end
indEnd = length(dTime);
if   sum([options.sources(:).type]==0) > 0
    if options.OnLine
        sigmaHat = posterior.a_sigma(:,dTime)./posterior.b_sigma(:,dTime);
        var_sigma = sigmaHat./posterior.b_sigma(dTime);
    else
        sigmaHat = posterior.a_sigma./posterior.b_sigma;
        var_sigma = sigmaHat./posterior.b_sigma;
    end
end
if options.dim.n > 0
    mux = posterior.muX(:,dTime);
    try
        vx = VBA_getVar(posterior.SigmaX.current,indEnd);
    catch
        vx = zeros(size(mux));
    end
    if ~any(isinf(posterior.a_alpha))
        if options.OnLine
            alphaHat = posterior.a_alpha(dTime)./posterior.b_alpha(dTime);
            var_alpha = alphaHat./posterior.b_alpha(dTime);
        else
            alphaHat = posterior.a_alpha./posterior.b_alpha;
            var_alpha = alphaHat./posterior.b_alpha;
        end
    else
        alphaHat = Inf;
        var_alpha = 0;
    end
    vx0 = VBA_getVar(posterior.SigmaX0);
    if options.updateX0
        dx0 = suffStat.dx0;
    else
        dx0 = posterior.muX0;
    end
end
if options.dim.n_theta > 0
    if options.OnLine
        dtheta = suffStat.dtheta(:,dTime);
    else
        dtheta = suffStat.dtheta;
    end
    vtheta = VBA_getVar(posterior.SigmaTheta,indEnd);
end
if options.dim.n_phi > 0
    if options.OnLine
        dphi = suffStat.dphi(:,dTime);
    else
        dphi = suffStat.dphi;
    end
    vphi = VBA_getVar(posterior.SigmaPhi,indEnd);
end

%% Vertical time encoding
if isequal(dTime,1) && size(y,1) > 1
    
    n_s = numel(options.sources);
    n_t = max(cellfun(@numel,{options.sources.out}));
    
    new_y = nan(n_s,n_t);
    new_gx = nan(n_s,n_t);
    new_vy = nan(n_s,n_t);
    new_isYout = ones(n_s,n_t);
    
    for si=1:n_s
        s_idx = options.sources(si).out;
        
        n_t_s = numel(options.sources(si).out);
        
        new_y(si, 1:n_t_s) = y(s_idx)';
        new_gx(si, 1:n_t_s) = gx(s_idx)';
        new_vy(si, 1:n_t_s) = vy(s_idx)';
        
        try 
            new_isYout(si, 1:n_t_s) = options.isYout(s_idx)';
        catch
            new_isYout(si, 1:n_t_s) = 0;
        end
      
    end
    
    y = new_y;
    vy = new_vy;
    gx = new_gx;
    options.isYout = new_isYout;
    
    if options.dim.n > 0 
        mux = mux';
        vx = vx';
    end
    
    dTime = 1:n_t;

        
    for si=1:n_s
        options.sources(si).out = si;
    end
end


switch flag % What piece of the model to display?
    
    case 'Y' % observation only
        
        % update top subplots
        update_observation_plot()
        
    case 'X' % Hidden-states related quantities
        if options.dim.n == 0, return; end

        % update top subplots
        update_observation_plot()
        
        % --
        
        % get display indices if delay embedding
        if sum(options.delays) > 0
            ind = 1:options.inF.dim.n;
        else
            ind = 1:size(mux,1);
        end
        
        % update middle-left subplot: hidden states
        %cla(display.ha(3))
        try
            plotUncertainTimeSeries(mux,vx,dTime,display.ha(3),ind);
        catch
            plotUncertainTimeSeries(mux,vx,[],display.ha(3),ind);
        end
        
        % update middle-right subplot: initial conditions
        if options.updateX0
            %cla(display.ha(4))
            plotUncertainTimeSeries(-dx0,vx0,1,display.ha(4));
        elseif isequal(it,0)
            plotUncertainTimeSeries(dx0,vx0,1,display.ha(4));
        end
        
        displayDF(F,display)
        
    case 'phi' % Observation parameters
        if options.dim.n_phi == 0, return; end

        % update top subplots
        update_observation_plot()
        
        % --
        
        % update bottom-left subplot: observation parameters
        if size(dphi,2) == 1 % for on-line wrapper
            dTime = 1;
        end
        %cla(display.ha(5))
        plotUncertainTimeSeries(-dphi,vphi,dTime,display.ha(5));
        
        displayDF(F,display)
        
    case 'theta' % Evolution parameters
        if options.dim.n_theta == 0, return; end
       
        % update bottom-right subplot: observation parameters
        if size(dtheta,2) == 1 % for on-line wrapper
            dTime = 1;
        end
        %cla(display.ha(7))
        plotUncertainTimeSeries(-dtheta,vtheta,dTime,display.ha(7));
        
        displayDF(F,display)   
        
    case 'precisions' % Precision hyperparameters
        
        % update top subplots
        update_observation_plot()
        
        % --
        
        if (options.updateHP || isequal(it,0)) && sum([options.sources(:).type]==0)>0
            dTime = 1;
            %cla(display.ha(6))
            logCI = (log(sigmaHat+sqrt(var_sigma)) - log(sigmaHat))';
            plotUncertainTimeSeries(log(sigmaHat'),logCI.^2,dTime,display.ha(6));
        end
        
        % update middle-right subplot: state noise
        if options.dim.n > 0 && ~any(isinf(alphaHat))
            dTime = 1;
            %cla(display.ha(8))
            logCI = log(alphaHat+sqrt(var_alpha)) - log(alphaHat);
            plotUncertainTimeSeries(log(alphaHat),logCI.^2,dTime,display.ha(8));
        end
        
        displayDF(F,display)
        
    case 'F' % Free energy
        
        % Output in main matlab window
        dF = diff(F);
        if it > 0 && options.verbose
            fprintf(['VB iteration #',...
                num2str(it),...
                '         F=','%e',...
                '         ... dF=','%4.3e'],F(end),dF(end))
            fprintf('\n')
        end
        
end

drawnow

% #########################################################################
% update_observation_plot ()
%
% update the top plots which shows the observations and model predictions 
% for the current source
% #########################################################################
function update_observation_plot ()
    
    % extract relevant values to be displayed
    % =====================================================================
    
    % get index of observations from current source
    s_out = options.sources(currentSource).out;
        
    % get observation for selected source
    y_src = y(s_out,:);
    y_src_in = y_src;
    y_src_in(options.isYout(s_out, :) == 1) = nan;
    y_src_out = y_src;
    y_src_out(options.isYout(s_out, :) == 0) = nan;
    
    % get predictions for selected source
    g_src = gx(s_out,:);
    g_src_in = g_src;
    g_src_in(options.isYout(s_out, :) == 1) = nan;
    g_src_out = g_src;
    g_src_out(options.isYout(s_out, :) == 0) = nan;
    
    % top left plot: timseries
    % =====================================================================
   
    % gaussian or binary source case: line + points
    % ---------------------------------------------------------------------
    if options.sources(currentSource).type < 2
            
        % predictive density
        vy_s = vy(s_out, :);
        resetColors (display.ha(1));
        plotUncertainTimeSeries (g_src(:, dTime), vy_s(:, dTime), dTime, display.ha(1));

        % data points
        p_in = findobj(display.ha(1),'Tag','yPoint');
        if ~ isempty(p_in)
            for i = 1 : numel (p_in)
                set(p_in(i), 'XData', dTime, 'YData', y_src_in(i,:));
            end
        else
            resetColors(display.ha(1));
            plot(display.ha(1),dTime,y_src_in','.','MarkerSize',9,'Tag','yPoint');
        end
            
        % excluded data points
        p_out = findobj(display.ha(1),'Tag','yOut');
        if ~ isempty(p_out)
            set(p_out,'XData',dTime, 'YData', y_src_out );
        else
            plot(display.ha(1),dTime, y_src_out, '.', 'MarkerEdgeColor',[.7 .7 .7],'MarkerSize',9,'Tag','yOut');
        end
            
        % data lines
        p_l = findobj(display.ha(1),'Tag','yLine');
        if ~ isempty(p_l)
            for i = 1 : numel (p_l)
                set(p_l(i), 'XData', dTime, 'YData', y_src(i,:));
            end
        else
            resetColors(display.ha(1));
            plot(display.ha(1),dTime,y_src',':','Tag','yLine','MarkerSize',9);
        end
                             
    % categorical source case: heatmap + points
    % ---------------------------------------------------------------------
    else
        
        % predictive density
        p_pred = findobj (display.ha(1), 'Tag', 'Ypred');
        if ~ isempty (p_pred)
            set (p_pred, 'CData', g_src);
        else
            imagesc (g_src, 'Parent', display.ha(1), 'Tag', 'yPred');
            set (display.ha(1), 'Clim', [0 1]) ;
            colormap (flipud (colormap ('bone')));   
        end
        
        % data points
        p_in = findobj (display.ha(1), 'Tag', 'yPoint');
        if ~ isempty (p_in)
            set (p_in, 'YData', VBA_indicator (y_src_in, [], true));
        else
            plot (display.ha(1), VBA_indicator (y_src_in, [], true), '.','ZData',2*ones(size(dTime)),'MarkerEdgeColor',[.9 0 0], 'MarkerSize', 9, 'Tag', 'yPoint');
        end
        
        % excluded points
        p_out = findobj (display.ha(1), 'Tag', 'yOut');
        if ~ isempty (p_in)
            set (p_out, 'YData', VBA_indicator (y_src_out, [], true));
        else
            plot (display.ha(1), VBA_indicator (y_src_out, [], true), '.','ZData',ones(size(dTime)),'MarkerEdgeColor',[.7 .7 .7], 'Tag', 'yOut');
        end
           
    end
    
    % ensure proper scaling of the axis
    set(display.ha(1),'XLim',[dTime(1)-0.5, dTime(end)+0.5]);
      
    % top right plot: predicted vs observed
    % =====================================================================
  
    % plot identity line
    % ---------------------------------------------------------------------
    % define boundaries
    if options.sources(currentSource).type == 0
    	miy = min (min ([g_src; y_src]));
        may = max (max ([g_src; y_src]));
    else
        miy = 0;
        may = 1;
    end
        
    % draw line
    hr = findobj (display.ha(2), 'Tag', 'refline');
    if ~ isempty (hr)
    	set (hr, 'XData',[miy, may], 'YData', [miy, may]);
    else
        plot (display.ha(2), [miy, may], [miy, may], 'k:', 'Tag', 'refline');
    end
   
    % gaussian source case
    % ---------------------------------------------------------------------
    if options.sources(currentSource).type == 0
     
        % included points
        p_in = findobj(display.ha(2), 'Tag', 'yIn');
        if ~ isempty(p_in)
            for i = 1 : numel (s_out)
                set (p_in(i), 'XData', g_src_in(i,:), 'YData', y_src_in(i,:));
            end
        else
            resetColors (display.ha(2));
            plot (display.ha(2), g_src_in', y_src_in', '.', 'MarkerSize', 9,'Tag', 'yIn');
        end

        % excluded points
        p_out = findobj (display.ha(2), 'Tag', 'yOut');
        if ~ isempty(p_out)
            set (p_out, 'XData', g_src_out(:), 'YData', y_src_out(:));
        else
            plot (display.ha(2), g_src_out(:), y_src_out(:), '.', 'MarkerEdgeColor', [.7 .7 .7], 'MarkerSize', 9, 'Tag', 'yOut');
        end
            
    % binary or categorical source case
    % ---------------------------------------------------------------------
    else
            
        % ellipse
        p_e = findobj(display.ha(2), 'Tag', 'ellipse');
        if isempty (p_e)
            grid_p = linspace (0, 1, 100);
            grid_std = sqrt (grid_p .* (1 - grid_p)); 
        
            plot (display.ha(2), grid_p, grid_p + grid_std, 'r--', 'Tag', 'ellpise');
            plot (display.ha(2), grid_p, grid_p - grid_std, 'r--', 'Tag', 'ellpise');
        end
        
        % included points
        [stackyin, stdyin, gridgin] = VBA_Bin2Cont (g_src_in, y_src_in);
        
        p_in = findobj (display.ha(2), 'Tag', 'yIn');
        if ~ isempty (p_in)
            set (p_in, ...
            'XData', gridgin, ...
            'YData', stackyin, ...
            'LData', stdyin, ...
            'UData', stdyin);
        else
            errorbar (gridgin, stackyin, stdyin, '.', 'MarkerSize', 9, 'parent', display.ha(2), 'Tag', 'yIn');
        end
        
        % excluded points
        [stackyout, stdyout, gridgout] = VBA_Bin2Cont (g_src_out, y_src_out);
         
        p_out = findobj (display.ha(2), 'Tag', 'yOut');
        if ~ isempty (p_out)
            set(p_out, ...
            'XData', gridgout, ...
            'YData', stackyout, ...
            'LData', stdyout, ...
            'UData', stdyout);
        else
            errorbar (gridgout, stackyout, stdyout, '.', 'Color', [.7 .7 .7], 'MarkerSize', 9, 'parent', display.ha(2), 'Tag', 'yOut');
        end   
    end       
end

% #########################################################################
% displayDF(F, display)
%
% update bottom text showing free energy updates
% #########################################################################
function displayDF(F,display)
    if ~ display.OnLine
        try
            dF = diff (F);
            set (display.ho, 'string', ['Model evidence: log p(y|m) >= ', num2str(F(end), '%1.3e'), ' , dF= ', num2str(dF(end), '%4.3e')])
        catch
            try
                set (display.ho, 'string', ['Model evidence: log p(y|m) >= ', num2str(F(end), '%4.3e')])
            end
        end
    end
end

% #########################################################################
% resetColors (h)
%
% restart color sequence before adding new data to a plot
% #########################################################################
function resetColors (h)
    if ~ verLessThan('matlab','8.4')
        set(h,'ColorOrderIndex',1);
    end
end

end
