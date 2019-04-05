function VBA_updateDisplay(posterior,suffStat,options,y,it,flag)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% VBA_updateDisplay(posterior,suffStat,options,y,it,flag)
% updates display of sufficient statistics
%
% This function deals with the screen display of iterative sufficient
% statistics updates of the VBA inversion algorithm

% IN:
%   - posterior: structure defining parmeter distribution
%   - suffStat: structure defining sufficient statistics
%   - options: usual options structure
%   - y: observations
%   - it: current iteration of the inversion
%   - flag: cell array of strings defining which plots to update.
%           by default, flag = {'X', 'phi', 'theta', 'precisions'} (update
%           all plots)
%
% /////////////////////////////////////////////////////////////////////////

% skip if no display requested
% =========================================================================
if ~ options.DisplayWin
    return
end

% dirty fix for online inversion
% =========================================================================
% preventing bar display on first iteration
if options.OnLine && it <= 1
    return
end
% index of last computed element
indEnd = size (y, 2);

% check flag
% =========================================================================
if ~ exist ('flag', 'var')
    % by default, update all plots
    flag = {'X', 'phi', 'theta', 'precisions'};
elseif ~ iscellstr (flag)
    flag = {flag};
end

% if this is a deterministic dynamical system inversion, expand the
% dynamics and display the result
% =========================================================================
if isequal (options.g_fname, @VBA_odeLim)  
    % rebuild posterior from dummy 'ODE' posterior
    [posterior, options, ~, suffStat] = VBA_odeLim2NLSS (posterior, options, options.dim, suffStat, [], 0);
    
    % call VBA_updateDisplay again
    VBA_updateDisplay(posterior,suffStat,options,y,it);
    
    % do not plot the compact system
    return
end

% get user actions
% =========================================================================

% check whether 'pause' button is toggled on
VBA_pause(options);

% replace source selector callback
ud = get(getPanel(options.display.hfp),'userdata') ;

    function update_plot ()
        cla(options.display.ha(1));
        cla(options.display.ha(2));
        VBA_updateDisplay(posterior,suffStat,options,y,0,{});
    end

ud.update_plot = @update_plot;
set(getPanel(options.display.hfp),'userdata',ud) ;

% get selected source
ud = VBA_check_struct(ud,'currentSource', 1 );
currentSource = ud.currentSource;

% =========================================================================
% update top subplots
update_observation_plot();

% Hidden-states related quantities
% =====================================================================
if ismember ('X', flag) && options.dim.n > 0
    
    % compute sufficient statistics to display
    mux = posterior.muX;
    try
        vx = VBA_getVar (posterior.SigmaX.current, indEnd);
    catch
        vx = zeros (size (mux));
    end
    % deal with vertical case
    if size (mux, 2) == 1 && size (mux, 1) > 1
            mux = mux';
            vx = vx';
    end
    % timeline
    T = size (mux, 2);
    % get display indices if delay embedding
    if sum(options.delays) > 0
        ind = 1 : options.inF.dim.n;
    else
        ind = 1 : size (mux, 1);
    end
    % update middle-left subplot: hidden states
    plotUncertainTimeSeries(mux,vx,1:T,options.display.ha(3),ind);
    % ensure proper scaling of the axis
    set(options.display.ha(3),'XLim',[0.5, T+0.5]);
    
    % update middle-right subplot: initial conditions
    vx0 = VBA_getVar (posterior.SigmaX0);
    if options.updateX0
        dx0 = suffStat.dx0;
    else
        dx0 = - posterior.muX0;
    end
    plotUncertainTimeSeries(- dx0,vx0,1,options.display.ha(4));
end

% Observation parameters
% =====================================================================
if ismember ('phi', flag) && options.dim.n_phi > 0
    
    % compute sufficient statistics to display
    if options.OnLine
        dphi = suffStat.dphi(:, indEnd);
    else
        dphi = suffStat.dphi;
    end
    vphi = VBA_getVar (posterior.SigmaPhi, indEnd);
    
    % update middle-left subplot: observation parameters
    plotUncertainTimeSeries(- dphi, vphi, 1, options.display.ha(5));

end

% Evolution parameters
% =====================================================================
if ismember ('theta', flag) && options.dim.n_theta > 0
    
    % compute sufficient statistics to display
    if options.OnLine
        dtheta = suffStat.dtheta(:, indEnd);
    else
        dtheta = suffStat.dtheta;
    end
    vtheta = VBA_getVar (posterior.SigmaTheta, indEnd);
    
    % update bottom-left subplot: evolution parameters
    plotUncertainTimeSeries(- dtheta,vtheta,1,options.display.ha(7));
end

% Precision hyperparameters
% =====================================================================
if ismember ('precisions', flag)
    
    % update middle-right subplot: observation noise
    if (options.updateHP || isequal(it,0)) && sum([options.sources(:).type]==0)>0
        % compute sufficient statistics
        if options.OnLine
            sigmaHat = posterior.a_sigma(:, indEnd) ./ posterior.b_sigma(:, indEnd);
            var_sigma = sigmaHat ./ posterior.b_sigma(indEnd);
        else
            sigmaHat = posterior.a_sigma ./ posterior.b_sigma;
            var_sigma = sigmaHat ./ posterior.b_sigma;
        end
        % display
        logCI = (log (sigmaHat + sqrt (var_sigma)) - log (sigmaHat))';
        plotUncertainTimeSeries (log (sigmaHat'), logCI.^2, 1, options.display.ha(6));
    end
    
    % update middle-right subplot: state noise
    if options.dim.n > 0 && ~ any (isinf (posterior.a_alpha))
        if options.OnLine
            alphaHat = posterior.a_alpha(indEnd) ./ posterior.b_alpha(indEnd);
            var_alpha = alphaHat ./ posterior.b_alpha(indEnd);
        else
            alphaHat = posterior.a_alpha ./ posterior.b_alpha;
            var_alpha = alphaHat ./ posterior.b_alpha;
        end
        logCI = log(alphaHat+sqrt(var_alpha)) - log(alphaHat);
        plotUncertainTimeSeries(log(alphaHat),logCI.^2,1,options.display.ha(8));
    end
     
end

% update free energy display
displayDF(real(suffStat.F),options.display)

drawnow

%% ########################################################################
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
        
        y_src = y(s_out,:);
        vy_src = suffStat.vy(s_out, :);
        g_src = suffStat.gx(s_out,:);

        % deal with vertical case
        if size (y, 2) == 1 && size (y, 1) > 1
            y_src = y_src';
            g_src = g_src';
            vy_src = vy_src';
        end
        
        % get observation for selected source
        y_src_in = y_src;
        y_src_in(options.isYout(s_out, :) == 1) = nan;
        y_src_out = y_src;
        y_src_out(options.isYout(s_out, :) == 0) = nan;
        
        % get predictions for selected source
        g_src_in = g_src;
        g_src_in(options.isYout(s_out, :) == 1) = nan;
        g_src_out = g_src;
        g_src_out(options.isYout(s_out, :) == 0) = nan;    
        
        % top left plot: timseries
        % =====================================================================
        
        % timeline
        T = size (g_src, 2);

        % gaussian or binary source case: line + points
        % ---------------------------------------------------------------------
        if options.sources(currentSource).type < 2
            
            % predictive density
            resetColors (options.display.ha(1));
            plotUncertainTimeSeries (g_src, vy_src, 1:T, options.display.ha(1));
            
            % axis
            if options.sources(currentSource).type == 0
                ylim(options.display.ha(1),'auto');
            else
                ylim(options.display.ha(1),[-0.2 1.2]);
            end
            
            % data lines
            p_l = findobj(options.display.ha(1),'Tag','yLine');
            if ~ isempty(p_l)
                for i = 1 : numel (p_l)
                    set(p_l(i), 'XData', 1:T, 'YData', y_src(i,:));
                end
            else
                resetColors(options.display.ha(1));
                plot(options.display.ha(1),1:T,y_src',':','Tag','yLine','MarkerSize',9);
            end  
            
            % excluded data points
            p_out = findobj(options.display.ha(1),'Tag','yOut');
            if ~ isempty(p_out)
                for i = 1 : numel (p_out)
                    set(p_out(i),'XData',1:T, 'YData', y_src_out(i,:) );
                end
            else
                plot(options.display.ha(1),1:T, y_src_out, '.', 'MarkerEdgeColor',[.7 .7 .7],'MarkerSize',9,'Tag','yOut');
            end
            
            % data points
            p_in = findobj(options.display.ha(1),'Tag','yPoint');
            if ~ isempty(p_in)
                for i = 1 : numel (p_in)
                    set(p_in(i), 'XData', 1:T, 'YData', y_src_in(i,:));
                end
            else
                resetColors(options.display.ha(1));
                plot(options.display.ha(1),1:T,y_src_in','.','MarkerSize',9,'Tag','yPoint');
            end
            
        % categorical source case: heatmap + points
        % ---------------------------------------------------------------------
        else
            
            % predictive density
            p_pred = findobj (options.display.ha(1), 'Tag', 'Ypred');
            if ~ isempty (p_pred)
                set (p_pred, 'CData', g_src);
            else
                imagesc (g_src, 'Parent', options.display.ha(1), 'Tag', 'yPred');
                set (options.display.ha(1), 'Clim', [0 1]) ;
                colormap (flipud (colormap ('bone')));
            end
            
            % excluded points
            p_out = findobj (options.display.ha(1), 'Tag', 'yOut');
            if ~ isempty (p_out)
                set (p_out, 'YData', VBA_indicator (y_src_out, [], true));
            else
                plot (options.display.ha(1), VBA_indicator (y_src_out, [], true), '.','ZData',ones(size(1:T)),'MarkerEdgeColor',[.7 .7 .7], 'Tag', 'yOut');
            end
            
            % data points
            p_in = findobj (options.display.ha(1), 'Tag', 'yPoint');
            if ~ isempty (p_in)
                set (p_in, 'YData', VBA_indicator (y_src_in, [], true));
            else
                plot (options.display.ha(1), VBA_indicator (y_src_in, [], true), '.','ZData',2*ones(size(1:T)),'MarkerEdgeColor',[.9 0 0], 'MarkerSize', 9, 'Tag', 'yPoint');
            end
    
        end
        
        % ensure proper scaling of the axis
        set(options.display.ha(1),'XLim',[0.5, T+0.5]);
        
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
        hr = findobj (options.display.ha(2), 'Tag', 'refline');
        if ~ isempty (hr)
            set (hr, 'XData',[miy, may], 'YData', [miy, may]);
        else
            plot (options.display.ha(2), [miy, may], [miy, may], 'k:', 'Tag', 'refline');
        end
        
        % gaussian source case
        % ---------------------------------------------------------------------
        if options.sources(currentSource).type == 0
             
            % excluded points
            p_out = findobj (options.display.ha(2), 'Tag', 'yOut');
            if ~ isempty(p_out)
                set (p_out, 'XData', g_src_out(:), 'YData', y_src_out(:));
            else
                plot (options.display.ha(2), g_src_out(:), y_src_out(:), '.', 'MarkerEdgeColor', [.7 .7 .7], 'MarkerSize', 9, 'Tag', 'yOut');
            end
            
            % included points
            p_in = findobj(options.display.ha(2), 'Tag', 'yIn');
            if ~ isempty(p_in)
                for i = 1 : numel (p_in)
                    set (p_in(i), 'XData', g_src_in(i,:), 'YData', y_src_in(i,:));
                end
            else
                resetColors (options.display.ha(2));
                plot (options.display.ha(2), g_src_in', y_src_in', '.', 'MarkerSize', 9,'Tag', 'yIn');
            end
            
            % binary or categorical source case
            % ---------------------------------------------------------------------
        else
            
            % ellipse
            p_e = findobj(options.display.ha(2), 'Tag', 'ellipse');
            if isempty (p_e)
                grid_p = linspace (0, 1, 100);
                grid_std = sqrt (grid_p .* (1 - grid_p));
                
                plot (options.display.ha(2), grid_p, grid_p + grid_std, 'r--', 'Tag', 'ellpise');
                plot (options.display.ha(2), grid_p, grid_p - grid_std, 'r--', 'Tag', 'ellpise');
            end
            
            % excluded points
            [stackyout, stdyout, gridgout] = VBA_Bin2Cont (g_src_out, y_src_out);
            
            p_out = findobj (options.display.ha(2), 'Tag', 'yOut');
            if ~ isempty (p_out)
                set(p_out, ...
                    'XData', gridgout, ...
                    'YData', stackyout, ...
                    'LData', stdyout, ...
                    'UData', stdyout);
            else
                errorbar (gridgout, stackyout, stdyout, '.', 'Color', [.7 .7 .7], 'MarkerSize', 9, 'parent', options.display.ha(2), 'Tag', 'yOut');
            end
            
            % included points
            [stackyin, stdyin, gridgin] = VBA_Bin2Cont (g_src_in, y_src_in);
            
            p_in = findobj (options.display.ha(2), 'Tag', 'yIn');
            if ~ isempty (p_in)
                set (p_in, ...
                    'XData', gridgin, ...
                    'YData', stackyin, ...
                    'LData', stdyin, ...
                    'UData', stdyin);
            else
                errorbar (gridgin, stackyin, stdyin, '.', 'MarkerSize', 9, 'parent', options.display.ha(2), 'Tag', 'yIn');
            end
        end
    end

% #########################################################################
% displayDF(F, options.display)
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
