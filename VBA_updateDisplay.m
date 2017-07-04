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

ud = check_struct(ud, ...
    'currentSource', 1 );

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
        try
            cla(display.ha(3))
            plotUncertainTimeSeries(mux,vx,dTime,display.ha(3),ind);
        catch
            cla(display.ha(3))
            plotUncertainTimeSeries(mux,vx,[],display.ha(3),ind);
        end
        set(display.ha(3),'ygrid','on','xgrid','off')
        axis(display.ha(3),'tight')
        
        % update middle-right subplot: initial conditions
        if options.updateX0
            cla(display.ha(4))
            plotUncertainTimeSeries(-dx0,vx0,1,display.ha(4));
            set(display.ha(4),'ygrid','on','xgrid','off')
        elseif isequal(it,0)
            plotUncertainTimeSeries(dx0,vx0,1,display.ha(4));
            set(display.ha(4),'ygrid','on','xgrid','off')
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
        set(display.ha(5),'ygrid','on','xgrid','off')
        
        displayDF(F,display)
        
    case 'theta' % Evolution parameters
        
        % update bottom-right subplot: observation parameters
        
        if size(dtheta,2) == 1 % for on-line wrapper
            dTime = 1;
        end
        cla(display.ha(7))
        plotUncertainTimeSeries(-dtheta,vtheta,dTime,display.ha(7));
        set(display.ha(7),'ygrid','on','xgrid','off')
        
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
            set(display.ha(6),'ygrid','on','xgrid','off')
        end
        
        % update middle-right subplot: state noise
        if options.dim.n > 0 && ~any(isinf(alphaHat))
            dTime = 1;
            cla(display.ha(8))
            logCI = log(alphaHat+sqrt(var_alpha)) - log(alphaHat);
            plotUncertainTimeSeries(log(alphaHat),logCI.^2,dTime,display.ha(8));
            set(display.ha(8),'ygrid','on','xgrid','off')
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
            plot(display.ha(1),dTime,y_s',':')
            plot(display.ha(1),dTime,y_s','.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',10)
            plot(display.ha(1),dTime,y_s_on','.','MarkerEdgeColor',[0 0 .8],'MarkerSize',10)
            vy_s= vy(s_out,:);
            plotUncertainTimeSeries(gx(s_out,dTime),vy_s(:,dTime),dTime,display.ha(1));
        else
            imagesc(gx(s_out,:),'Parent',display.ha(1));
            set(display.ha(1),'Clim',[0 1]) ;
            colormap(flipud(colormap('bone')));
            plot(display.ha(1),multi2num(y_s_on),'.r');
        end
        set(display.ha(1),'ygrid','on','xgrid','off')
        axis(display.ha(1),'tight')
        
        % update top-right subplot: predicted VS observed data
        cla(display.ha(2))
        % plot identity
        if options.sources(currentSource).type==0
            miy = min(min([gx(s_out,:);y(s_out,:)]));
            may = max(max([gx(s_out,:);y(s_out,:)]));
            plot(display.ha(2),[miy,may],[miy,may],'r')
        else
            plot(display.ha(2),[0,1],[0,1],'r')
        end
        
        if options.sources(currentSource).type==0
            gx_src = gx(s_out,:) ;
            y_src = y(s_out,:) ;
            gxout = gx_src(~~options.isYout(s_out,:));
            yout = y_src(~~options.isYout(s_out,:));
            gxin = gx_src(~options.isYout(s_out,:));
            yin = y_src(~options.isYout(s_out,:));
            plot(display.ha(2),gxout(:),yout(:),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',10)
            plot(display.ha(2),gxin(:),yin(:),'.','MarkerEdgeColor',[0 0 .8],'MarkerSize',10)
        else
            gridp = 0:1e-2:1;
            plot(display.ha(2),gridp,gridp+sqrt(gridp.*(1-gridp)),'r--')
            plot(display.ha(2),gridp,gridp-sqrt(gridp.*(1-gridp)),'r--')
            errorbar(gridgout{currentSource},stackyout{currentSource},stdyout{currentSource},'.','Color',[.85 .85 .85],'parent',display.ha(2))
            errorbar(gridgin{currentSource},stackyin{currentSource},stdyin{currentSource},'k.','parent',display.ha(2))
        end
        
        grid(display.ha(2),'on')
        axis(display.ha(2),'tight')
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

    function [nseq] = interpOut(seq,isYout)
        nseq=seq;
        return
        [n,t]=size(seq);
        times=1:t;
        for i = 1:n
            mask =  isYout(i,:)==0;
            nseq_i = seq(i,:);
            nseq_i(~mask) = interp1(times(mask), nseq_i(mask), times(~mask),'spline',0);
            nseq_i(isnan(nseq_i))=0;
            nseq(i,:)=nseq_i;
        end
        
    end

    function updateSource(hObject,evt,si)
        changeSource(hObject,evt,si);
        VBA_updateDisplay(posterior,suffStat,options,y,0,'Y');
    end

end

