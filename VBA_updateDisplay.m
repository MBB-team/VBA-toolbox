function VBA_updateDisplay(posterior,suffStat,options,y,it,flag)
% updates display of sufficient statistics
% function VBA_updateDisplay(F,posterior,suffStat,options,y,it,display,flag)
% This function deals with the screen display of iterative sufficient
% statistics updates of the VBA inversion algorithm

if ~options.DisplayWin
    return
end

if true %numel(options.sources)>1 || options.sources(1).type==2
    VBA_updateDisplay_extended(posterior,suffStat,options,y,it,flag);
    %fprintf('Display extended\n');
    return
end

F = real(suffStat.F); % should not be required (cf. exit sanity check in VBA_NLStateSpaceModel.m)!
display = options.display;

% check whether 'pause' button is toggled on
VBA_pause(options)

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
dTime = [1:size(y,2)];
try
    gx = suffStat.gx(:,dTime);
    vy = suffStat.vy(:,dTime);
catch
    gx = [];
    vy = [];
end
indEnd = length(dTime);
if  ~options.binomial
    if options.OnLine
        sigmaHat = posterior.a_sigma(dTime)./posterior.b_sigma(dTime);
        var_sigma = sigmaHat./posterior.b_sigma(dTime);
    else
        sigmaHat = posterior.a_sigma./posterior.b_sigma;
        var_sigma = sigmaHat./posterior.b_sigma;
    end
else
    try
        [stackyin,stdyin,gridgin] = VBA_Bin2Cont(gx(~options.isYout),y(~options.isYout));
        [stackyout,stdyout,gridgout] = VBA_Bin2Cont(gx(~~options.isYout),y(~~options.isYout));
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
    dTime = [1:size(y,2)];
end



switch flag % What piece of the model to display?
   
   
    case 'X' % Hidden-states related quantities
       
       
        % update top-left subplot: predictive density
        cla(display.ha(1))
        if ~isempty(gx) && ~isempty(vy)
            plot(display.ha(1),dTime,y','LineStyle',':','marker','.')
            plotUncertainTimeSeries(gx,vy,dTime,display.ha(1));
            set(display.ha(1),'ygrid','on','xgrid','off')
            axis(display.ha(1),'tight')
        end
        
        % update top-right subplot: predicted VS observed data
        cla(display.ha(2))
        if ~isempty(gx) && ~isempty(vy)
            if  ~options.binomial
                miy = min([gx(:);y(:)]);
                may = max([gx(:);y(:)]);
                plot(display.ha(2),[miy,may],[miy,may],'r')
                gxout = gx(~~options.isYout);
                yout = y(~~options.isYout);
                gxin = gx(~options.isYout);
                yin = y(~options.isYout);
                plot(display.ha(2),gxout(:),yout(:),'r.')
                plot(display.ha(2),gxin(:),yin(:),'k.')
                if ~isempty(yout)
                    legend(display.ha(2),{'','excluded','fitted'})
                end
            else
                plot(display.ha(2),[0,1],[0,1],'r')
                gridp = 0:1e-2:1;
                plot(display.ha(2),gridp,gridp+sqrt(gridp.*(1-gridp)),'r--')
                plot(display.ha(2),gridp,gridp-sqrt(gridp.*(1-gridp)),'r--')
                errorbar(gridgout,stackyout,stdyout,'r.','parent',display.ha(2))
                errorbar(gridgin,stackyin,stdyin,'k.','parent',display.ha(2))
                if ~isempty(gridgout)
                    legend(display.ha(2),{'','','','excluded','fitted'})
                end
            end
            grid(display.ha(2),'on')
            axis(display.ha(2),'tight')
        end
        
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
        
        
        % update top-left subplot: predictive density
        cla(display.ha(1))
        plot(display.ha(1),dTime,y',':')
        plot(display.ha(1),dTime,y','.')
        if ~isempty(gx) && ~isempty(vy)
            plotUncertainTimeSeries(gx,vy,dTime,display.ha(1));
        end
        set(display.ha(1),'ygrid','on','xgrid','off')
        axis(display.ha(1),'tight')
        
        % update top-right subplot: predicted VS observed data
        cla(display.ha(2))
        if ~isempty(gx) && ~isempty(vy)
            if  ~options.binomial
                miy = min([gx(:);y(:)]);
                may = max([gx(:);y(:)]);
                plot(display.ha(2),[miy,may],[miy,may],'r')
                gxout = gx(~~options.isYout);
                yout = y(~~options.isYout);
                gxin = gx(~options.isYout);
                yin = y(~options.isYout);
                plot(display.ha(2),gxout(:),yout(:),'r.')
                plot(display.ha(2),gxin(:),yin(:),'k.')
                if ~isempty(yout)
                    legend(display.ha(2),{'','excluded','fitted'})
                end
            else
                plot(display.ha(2),[0,1],[0,1],'r')
                gridp = 0:1e-2:1;
                plot(display.ha(2),gridp,gridp+sqrt(gridp.*(1-gridp)),'r--')
                plot(display.ha(2),gridp,gridp-sqrt(gridp.*(1-gridp)),'r--')
                errorbar(gridgout,stackyout,stdyout,'r.','parent',display.ha(2))
                errorbar(gridgin,stackyin,stdyin,'k.','parent',display.ha(2))
                if ~isempty(gridgout)
                    legend(display.ha(2),{'','','','excluded','fitted'})
                end
            end
            grid(display.ha(2),'on')
            axis(display.ha(2),'tight')
        end
        
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
        
        % update top-left subplot: predictive density
        cla(display.ha(1))
        plot(display.ha(1),dTime,y',':')
        plot(display.ha(1),dTime,y','.')
        if ~isempty(gx) && ~isempty(vy)
            plotUncertainTimeSeries(gx,vy,dTime,display.ha(1));
        end
        
        
        if options.updateHP || (isequal(it,0) && ~options.binomial)
            
            % update middle-left subplot: measurement noise
            if ~options.binomial
                if size(sigmaHat,2) > 1  % for on-line wrapper
                    cla(display.ha(6))
                else
                    dTime = it+1;
                    set(display.ha(6),'xlim',[.2,it+1.8],'xtick',[])
                end
                logCI = log(sigmaHat+sqrt(var_sigma)) - log(sigmaHat);
                plotUncertainTimeSeries(log(sigmaHat),logCI.^2,dTime,display.ha(6));
                set(display.ha(6),'ygrid','on','xgrid','off')
            end
            
            % update middle-right subplot: state noise
            if options.dim.n > 0 && ~any(isinf(alphaHat))
                if size(alphaHat,2) > 1  % for on-line wrapper
                    cla(display.ha(8))
                else
                    dTime = it+1;
                    set(display.ha(8),'xlim',[.2,it+1.8],'xtick',[])
                end
                logCI = log(alphaHat+sqrt(var_alpha)) - log(alphaHat);
                plotUncertainTimeSeries(log(alphaHat),logCI.^2,dTime,display.ha(8));
                set(display.ha(8),'ygrid','on','xgrid','off')
            end
            
            displayDF(F,display)
            
        end
        
    case 'F' % Free energy
        
        % Output in main matlab window
        dF = diff(F);
        if it > 0 && options.verbose
            fprintf(['VB iteration #',num2str(it),'         F=','%e','         ... dF=','%4.3e'],F(end),dF(end))
            fprintf('\n')
        end
        
end

%drawnow


%--- subfunction ---%

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

