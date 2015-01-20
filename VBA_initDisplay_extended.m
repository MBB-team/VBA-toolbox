function [options] = VBA_initDisplay_extended(options,priors)
% initialize the graphical output

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

try % if called from uitabs of VBA_ReDisplay
    display.hfp = options.display.hfp;
    get(display.hfp,'children');
catch
    display.hfp = findobj('tag','VBNLSS');
end
if isempty(display.hfp)
    pos0 = get(0,'screenSize');
    pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];
    display.hfp = figure('position',pos,'color',[1 1 1],'name',options.figName,'menubar','none','tag','VBNLSS','Renderer','OpenGL','visible',visible);
else
    display.hfp = display.hfp(1);
    hc = intersect(findobj('tag','VBLaplace'),get(display.hfp,'children'));
    if ~isempty(hc)
        delete(hc)
    end
    set(display.hfp,'name',options.figName,'visible',visible);
end

Ns = numel(options.sources);


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

hPanel = getPanel(display.hfp);

for s_i=1:Ns
    
display.ha(2*s_i-1) = subplot(3+Ns,2,2*s_i-1,'parent',hPanel,'xlim',xlim,'nextplot','add','tag','VBLaplace','box','off');
if ~priors
    title(display.ha(2*s_i-1),'posterior predictive density: p(g(x)|y,m)','fontsize',11)
else
    title(display.ha(2*s_i-1),'prior predictive density: p(g(x)|m)','fontsize',11)
end
xlabel(display.ha(2*s_i-1),xl,'fontsize',8)
if ~priors
    ylabel(display.ha(2*s_i-1),'<g(x)|y,m> & y','fontsize',8)
else
    ylabel(display.ha(2*s_i-1),'<g(x)|m> & y','fontsize',8)
end

display.ha(2*s_i) = subplot(3+Ns,2,2*s_i,'parent',hPanel,'nextplot','add','tag','VBLaplace','box','off');
if ~priors
    title(display.ha(2*s_i),'Model fit: <g(x)|y,m> versus y','fontsize',11)
    xlabel(display.ha(2*s_i),'<g(x)|y,m>','fontsize',8)
else
    title(display.ha(2*s_i),'Model fit: <g(x)|m> versus y','fontsize',11)
    xlabel(display.ha(2*s_i),'<g(x)|m>','fontsize',8)
end
ylabel(display.ha(2*s_i),'y','fontsize',8)

end
% Create axes for hidden states and initial conditions
if options.dim.n > 0
    display.ha(2*Ns+1) = subplot(3+Ns,2,2*Ns+1,'parent',hPanel,'nextplot','add','tag','VBLaplace','box','off');
    if ~priors
        title(display.ha(2*Ns+1),'hidden states: p(x|y,m)','fontsize',11)
    else
        title(display.ha(2*Ns+1),'hidden states: p(x|m)','fontsize',11)
    end
    xlabel(display.ha(2*Ns+1),'time','fontsize',8)
    if ~priors
        ylabel(display.ha(2*Ns+1),'<x|y,m>','fontsize',8)
    else
        ylabel(display.ha(2*Ns+1),'<x|m>','fontsize',8)
    end
    display.ha(2*Ns+2) = subplot(3+Ns,2,2*Ns+2,'parent',hPanel,'nextplot','add','xlim',[0.2,options.dim.n+0.8],'xtick',[],'tag','VBLaplace','box','off');
    if ~priors
        title( display.ha(2*Ns+2),'initial conditions: p(x_0|y,m)','fontsize',11)
    else
        title( display.ha(2*Ns+2),'initial conditions: p(x_0|m)','fontsize',11)
    end
    if options.updateX0
        xlabel(display.ha(2*Ns+2),'x_0 dimensions','fontsize',8)
        if ~priors
            ylabel(display.ha(2*Ns+2),'<x_0|y,m> - <x_0|m>','fontsize',8)
        else
            ylabel(display.ha(2*Ns+2),'<x_0|m>','fontsize',8)
        end
    else
        xlabel(display.ha(2*Ns+2),'x_0 dimensions [fixed pdf]','fontsize',8)
        ylabel(display.ha(2*Ns+2),'x_0','fontsize',8)
    end
end

% Create axes for observation parameters
if options.dim.n_phi > 0
    display.ha(2*Ns+3) = subplot(3+Ns,2,2*Ns+3,'parent',hPanel,'nextplot','add','xlim',[0.2,options.dim.n_phi+0.8],'xtick',[],'tag','VBLaplace','box','off');
    if ~priors
        title(display.ha(2*Ns+3),'observation parameters: p(phi|y,m)','fontsize',11)
    else
        title(display.ha(2*Ns+3),'observation parameters: p(phi|m)','fontsize',11)
    end
    if ~options.OnLine
        xlabel(display.ha(2*Ns+3),'phi dimensions','fontsize',8)
    else
        xlabel(display.ha(2*Ns+3),'time','fontsize',8)
    end
    if ~priors
        ylabel(display.ha(2*Ns+3),'<phi|y,m> - <phi|m>','fontsize',8)
    else
        ylabel(display.ha(2*Ns+3),'<phi|m>','fontsize',8)
    end
end

% Create axes for measurement noise precision hyperparameter
Ngs=sum([options.sources(:).type]==0);
if Ngs>0
    display.ha(2*Ns+4) = subplot(3+Ns,2,2*Ns+4,'parent',hPanel,'nextplot','add','xlim',[0.2,Ngs+0.8],'xtick',[],'tag','VBLaplace','box','off');
    if ~priors
        title(display.ha(2*Ns+4),'measurement noise precision: p(sigma|y,m)','fontsize',11)
    else
        title(display.ha(2*Ns+4),'measurement noise precision: p(sigma|m)','fontsize',11)
    end
    if ~options.OnLine && options.updateHP
        xlabel(display.ha(2*Ns+4),'channel','fontsize',8)
    elseif ~options.OnLine && ~options.updateHP
        xlabel(display.ha(2*Ns+4),'[fixed pdf]','fontsize',8)
    else
        xlabel(display.ha(2*Ns+4),'time','fontsize',8)
    end
    ylabel(display.ha(2*Ns+4),'<log(sigma)>','fontsize',8)
end

% Create axes for evolution parameters
if options.dim.n_theta > 0
    display.ha(2*Ns+5) = subplot(3+Ns,2,2*Ns+5,'parent',hPanel,'nextplot','add','xlim',[0.2,options.dim.n_theta+0.8],'xtick',[],'tag','VBLaplace','box','off');
    if ~priors
        title(display.ha(2*Ns+5),'evolution parameters: p(theta|y,m)','fontsize',11)
    else
        title(display.ha(2*Ns+5),'evolution parameters: p(theta|m)','fontsize',11)
    end
    if ~options.OnLine
        xlabel(display.ha(2*Ns+5),'theta dimensions','fontsize',8)
    else
        xlabel(display.ha(2*Ns+5),'time','fontsize',8)
    end
    if ~priors
        ylabel(display.ha(2*Ns+5),'<theta|y,m> - <theta|m>','fontsize',8)
    else
        ylabel(display.ha(2*Ns+5),'<theta|m>','fontsize',8)
    end
end

% Create axes for state noise precision hyperparameter
if ~isequal(options0.g_fname,@VBA_odeLim) && options.dim.n > 0 % not for non stochastic systems
    display.ha(2*Ns+6) = subplot(3+Ns,2,2*Ns+6,'parent',hPanel,'xlim',[0.2,1.8],'nextplot','add','tag','VBLaplace','box','off');
    if ~priors
        title(display.ha(2*Ns+6),'system''s noise precision: p(alpha|y,m)','fontsize',11)
    else
        title(display.ha(2*Ns+6),'system''s noise precision: p(alpha|m)','fontsize',11)
    end
    if ~options.OnLine && options.updateHP
        xlabel(display.ha(2*Ns+6),'VB iterations','fontsize',8)
    elseif ~options.OnLine && ~options.updateHP
        xlabel(display.ha(2*Ns+6),'[fixed pdf]','fontsize',8)
    else
        xlabel(display.ha(2*Ns+6),'time','fontsize',8)
    end
    ylabel(display.ha(2*Ns+6),'<log(alpha)>','fontsize',8)
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

drawnow
display.OnLine = options.OnLine;
options = options0;
options.display = display;

try
    getSubplots
end





