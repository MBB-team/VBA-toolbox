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
ud.update_plot = @() VBA_updateDisplay(posterior,suffStat,options,y,0,'Y');
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
    if dim.n_phi > 0
        VBA_updateDisplay(posterior,suffStat,options,y,it,'phi')
    end
    if dim.n > 0
        VBA_updateDisplay(posterior,suffStat,options,y,it,'X')
    end
    if dim.n_theta > 0
        VBA_updateDisplay(posterior,suffStat,options,y,it,'theta')
    end
    
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

% check time dimension
if isequal(dTime,1) && size(y,1) > 1
    gx = gx';
    vy = vy';
    y = y';
    if options.dim.n > 0
        mux = mux';
        vx = vx';
    end
    try
        options.isYout = options.isYout';
    end
    dTime = 1:size(y,2);
    for si=1:numel(options.sources)
        options.sources(si).out = si;
    end
end

Ns = numel(options.sources);
for s_i=1:Ns
    if  options.sources(s_i).type>0
        s_out = options.sources(s_i).out ;
        gx_out = gx(s_out,:);
        y_out = y(s_out,:);
        [stackyin{s_i},stdyin{s_i},gridgin{s_i}] = VBA_Bin2Cont(gx_out(~options.isYout(s_out,:)),y_out(~options.isYout(s_out,:)));
        [stackyout{s_i},stdyout{s_i},gridgout{s_i}] = VBA_Bin2Cont(gx_out(~~options.isYout(s_out,:)),y_out(~~options.isYout(s_out,:)));
    end
end

switch flag % What piece of the model to display?
    
    case 'Y' % observation only
        
        % update top subplots
        update_observation_plot()
        
    case 'X' % Hidden-states related quantities
        
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
        cla(display.ha(3))
        try
            plotUncertainTimeSeries(mux,vx,dTime,display.ha(3),ind);
        catch
            plotUncertainTimeSeries(mux,vx,[],display.ha(3),ind);
        end
        
        % update middle-right subplot: initial conditions
        if options.updateX0
            cla(display.ha(4))
            plotUncertainTimeSeries(-dx0,vx0,1,display.ha(4));
        elseif isequal(it,0)
            plotUncertainTimeSeries(dx0,vx0,1,display.ha(4));
        end
        
        displayDF(F,display)
        
    case 'phi' % Observation parameters
        
        % update top subplots
        update_observation_plot()
        
        % --
        
        % update bottom-left subplot: observation parameters
        if size(dphi,2) == 1 % for on-line wrapper
            dTime = 1;
        end
        cla(display.ha(5))
        plotUncertainTimeSeries(-dphi,vphi,dTime,display.ha(5));
        
        displayDF(F,display)
        
    case 'theta' % Evolution parameters
        
        % update bottom-right subplot: observation parameters
        if size(dtheta,2) == 1 % for on-line wrapper
            dTime = 1;
        end
        cla(display.ha(7))
        plotUncertainTimeSeries(-dtheta,vtheta,dTime,display.ha(7));
        
        displayDF(F,display)   
        
    case 'precisions' % Precision hyperparameters
        
        % update top subplots
        update_observation_plot()
        
        % --
        
        if (options.updateHP || isequal(it,0)) && sum([options.sources(:).type]==0)>0
            dTime = 1;
            cla(display.ha(6))
            logCI = (log(sigmaHat+sqrt(var_sigma)) - log(sigmaHat))';
            plotUncertainTimeSeries(log(sigmaHat'),logCI.^2,dTime,display.ha(6));
        end
        
        % update middle-right subplot: state noise
        if options.dim.n > 0 && ~any(isinf(alphaHat))
            dTime = 1;
            cla(display.ha(8))
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


%--- subfunction ---%

    function update_observation_plot()
        
        s_out=options.sources(currentSource).out;
        
        % update top-left subplot: predictive density
        cla(display.ha(1))
        y_s = y(s_out,:);
        y_s_on = y_s;
        y_s_on(options.isYout(s_out,:)==1)=nan;
        if options.sources(currentSource).type < 2
            p_l  = plot(display.ha(1),dTime,y_s',':');
            plot(display.ha(1),dTime,y_s','.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',9);
            p_mi = plot(display.ha(1),dTime,y_s_on','.','MarkerSize',9);
            
            vy_s= vy(s_out,:);
            [~,p_vr,p_vl] = plotUncertainTimeSeries(gx(s_out,dTime),vy_s(:,dTime),dTime,display.ha(1));
            for i=1:numel(p_l)
                set(p_mi(i),'MarkerEdgeColor',get(p_l(i),'Color'))
                set(p_vl(i),'Color',get(p_l(i),'Color'))
                set(p_vr(i),'FaceColor',get(p_l(i),'Color'))
            end
        else
            imagesc(gx(s_out,:),'Parent',display.ha(1));
            set(display.ha(1),'Clim',[0 1]) ;
            colormap(flipud(colormap('bone')));
            plot(display.ha(1),multi2num(y_s_on),'.r');
        end
        
        % update top-right subplot: predicted VS observed data
        cla(display.ha(2))
        % plot identity
        if options.sources(currentSource).type==0
            miy = min(min([gx(s_out,:);y(s_out,:)]));
            may = max(max([gx(s_out,:);y(s_out,:)]));
            plot(display.ha(2),[miy,may],[miy,may],'k:')
        else
            plot(display.ha(2),[0,1],[0,1],'k:')
        end
        
        if options.sources(currentSource).type==0
            gx_src = gx(s_out,:) ;
            y_src = y(s_out,:) ;
            gxout = gx_src(~~options.isYout(s_out,:));
            yout = y_src(~~options.isYout(s_out,:));
            gxin = gx_src;
            gxin(~~options.isYout(s_out,:)) = nan;
            yin = y_src;
            yin(~~options.isYout(s_out,:)) = nan;
            plot(display.ha(2),gxout(:),yout(:),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',9)
            plot(display.ha(2),gxin',yin','.','MarkerSize',9)
        else
            gridp = 0:1e-2:1;
            plot(display.ha(2),gridp,gridp+sqrt(gridp.*(1-gridp)),'r--')
            plot(display.ha(2),gridp,gridp-sqrt(gridp.*(1-gridp)),'r--')
            errorbar(gridgout{currentSource},stackyout{currentSource},stdyout{currentSource},'.','Color',[.7 .7 .7],'MarkerSize',9,'parent',display.ha(2))
            errorbar(gridgin{currentSource},stackyin{currentSource},stdyin{currentSource},'.','MarkerSize',9,'parent',display.ha(2))
        end
         
    end

    function [] = displayDF(F,display)
        if ~display.OnLine
            try
                dF = diff(F);
                set(display.ho,'string',['Model evidence: log p(y|m) >= ',num2str(F(end),'%1.3e'),' , dF= ',num2str(dF(end),'%4.3e')])
            catch
                try
                    set(display.ho,'string',['Model evidence: log p(y|m) >= ',num2str(F(end),'%4.3e')])
                end
            end
        end
    end

end

