function [hfp] = VBA_ReDisplay(posterior,out,newFig,fromPause)
% re-creates the graphical output of the VBA inversion + diagnostics
% function [hfp] = VBA_ReDisplay(posterior,out)
% NB: hfp is the figure handle. Note that VBA_ReDisplay first looks for a
% figure with a tag 'VBNLSS', i.e. a figure that was already opened to
% review a model inversion, and clears it if it finds it.


try; newFig; catch; newFig = 0; end
try; fromPause; catch; fromPause = 0; end

options = out.options;
options.DisplayWin = 1;
hfp = findobj('tag','VBNLSS');
if isempty(hfp) || newFig
    pos0 = get(0,'screenSize');
    pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];
    hfp = figure('position',pos,'color',[1 1 1],'name','VB-Laplace approximate Bayesian inference','menubar','none','tag','VBNLSS','Renderer','OpenGL');
else
    hfp = hfp(1);
    clf(hfp)
    set(hfp,'name',options.figName);
end

ud.posterior = posterior;
ud.out = out;
ud.diagnostics = getDiagnostics(posterior,out);

set(hfp,'userdata',ud);
if ~isempty(ud.diagnostics.kernels)
    if ~isinf(out.options.priors.a_alpha) && ~isequal(out.options.priors.b_alpha,0) && ~out.options.OnLine
        labels = {'summary','VB inversion','diagnostics','kernels','conv','priors','deterministic'};
        callbacks = {@mySummary,@myVB,@myDiagnostics,@myKernels,@myConv,@myPriors,@myDeterministic};
    else
        labels = {'summary','VB inversion','diagnostics','kernels','conv','priors'};
        callbacks = {@mySummary,@myVB,@myDiagnostics,@myKernels,@myConv,@myPriors};
    end
else
    labels = {'summary','VB inversion','diagnostics','conv','priors'};
    callbacks = {@mySummary,@myVB,@myDiagnostics,@myConv,@myPriors};
end
if fromPause
    active = 2;
else
    active = 1;
end
[handles] = spm_uitab(hfp,labels,callbacks,'diagnostics_tabs',active);
set(handles.htab,'backgroundcolor',[1 1 1])
set(handles.hh,'backgroundcolor',[1 1 1])
set(handles.hp,'HighlightColor',0.8*[1 1 1])
set(handles.hp,'backgroundcolor',[1 1 1])
set(handles.htab(1),'tooltipstring','summary description of the VB inversion')
set(handles.htab(2),'tooltipstring','results of the VB inversion (posterior pdfs)')
set(handles.htab(3),'tooltipstring','VB inversion diagnostics (residuals and parameters covariance matrices)')
if ~isempty(ud.diagnostics.kernels)
    if ~isinf(out.options.priors.a_alpha) && ~isequal(out.options.priors.b_alpha,0) && ~out.options.OnLine
        set(handles.htab(4),'tooltipstring','impulse response of the system to deterministic inputs')
        set(handles.htab(5),'tooltipstring','history of free energy values along VB optimization')
        set(handles.htab(6),'tooltipstring','priors and associated predictive densities (under the Laplace assumption)')
        set(handles.htab(7),'tooltipstring','results of the VB inversion of the deterministic system')
    else
        set(handles.htab(4),'tooltipstring','impulse response of the system to deterministic inputs')
        set(handles.htab(5),'tooltipstring','history of free energy values along VB optimization')
        set(handles.htab(6),'tooltipstring','priors and associated predictive densities (under the Laplace assumption)')
    end
else
    set(handles.htab(4),'tooltipstring','history of free energy values along VB optimization')
    set(handles.htab(5),'tooltipstring','priors and associated predictive densities (under the Laplace assumption)')
end
if fromPause
    feval(@myVB,hfp)
else
    feval(@mySummary,hfp)
end

function mySummary(hfp)

try
    hf = hfp;
catch
    hf = get(gco,'parent');
end
hc = intersect(findobj('tag','VBLaplace'),get(hf,'children'));
if ~isempty(hc)
    delete(hc)
end

ud = get(hf,'userdata');
out = ud.out;
diagnostics = ud.diagnostics;

try F = out.F(end); catch, F = '?'; end

str{1} = sprintf(['Date: ',datestr(out.date),'\n ']);
if ~out.options.OnLine
    s0 = ['VB converged in ',num2str(out.it),' iterations'];
else
    s0 = ['Online VB algorithm'];
end
try
    if floor(out.dt./60) == 0
        timeString = [num2str(floor(out.dt)),' sec'];
    else
        timeString = [num2str(floor(out.dt./60)),' min'];
    end
    str{2} = sprintf([s0,' (took ~',timeString,')','\n']);
catch
    str{2} = sprintf([s0,'\n']);
end
str{3} = sprintf(['Dimensions of the model:','\n ',...
    '    - data: p=',num2str(out.dim.p),'\n ',...
    '    - time samples: t=',num2str(out.dim.n_t),'\n ',...
    '    - hidden states: n=',num2str(out.dim.n),'\n ',...
    '    - evolution parameters: n_theta=',num2str(out.dim.n_theta),'\n ',...
    '    - observation parameters: n_phi=',num2str(out.dim.n_phi),'\n ']);
if out.options.binomial
    tmp = ' (binomial data)';
else
    tmp = [];
end
if out.dim.n >= 1
    if isinf(out.options.priors.a_alpha) && isequal(out.options.priors.b_alpha,0)
        str{4} = 'This was a deterministic dynamical system';
    else
        str{4} = 'This was a stochastic dynamical system';
    end
    if isa(out.options.g_fname,'function_handle')
        gfn = func2str(out.options.g_fname);
    else
        gfn = out.options.g_fname;
    end
    if isequal(gfn,'g_embed')
        gfn0 = out.options.inG.g_fname;
        if isa(gfn0,'function_handle')
            gfn0 = func2str(gfn0);
        end
        gfn = [gfn,' (',gfn0,')'];
        str{4} = [str{4},' (with delay embedding)'];
    end
    if isa(out.options.f_fname,'function_handle')
        ffn = func2str(out.options.f_fname);
    else
        ffn = out.options.f_fname;
    end
    if isequal(ffn,'f_embed')
        ffn0 = out.options.inF.f_fname;
        if isa(ffn0,'function_handle')
            ffn0 = func2str(ffn0);
        end
        ffn = [ffn,' (',ffn0,')'];
    end
    str{5} = sprintf(['    - observation function: ',gfn,tmp,'\n','    - evolution function: ',ffn,'\n ']);
else
    str{4} = 'The model was static (no hidden states)';
    if isa(out.options.g_fname,'function_handle')
        gfn = func2str(out.options.g_fname);
    else
        gfn = out.options.g_fname;
    end
    str{5} = sprintf(['    - observation function: ',gfn,tmp,'\n ']);
end
str{6} = ['Log model evidences:'];
str{7} = ['    - full model: log p(y|m) > ',num2str(F,'%4.3e')];
str{8} = ['    - null hypothesis: log p(y|H0) = ',...
    num2str(diagnostics.LLH0,'%4.3e')];
if ~out.options.OnLine && out.dim.n >= 1 && ~isinf(out.options.priors.a_alpha) && ~isequal(out.options.priors.b_alpha,0)
    Fd = out.options.init.out.F;
    str{9} = sprintf(['    - deterministic variant: log p(y|m,eta=0) > ',num2str(Fd,'%4.3e'),'\n ']);
else
    str{9} = [' '];
end
str{10} = sprintf(['Estimation efficiency (minus posterior entropies):','\n ']);
% str{11} = sprintf(['Information gain (Kullback-Leibler divergences DKL{prior||posterior}):','\n ']);
if ~isnan(diagnostics.efficiency.X)
    str{10} = sprintf([str{10},'    - hidden states: ',num2str(diagnostics.efficiency.X,'%4.3e'),'\n ']);
%     str{11} = sprintf([str{11},'    - hidden states: ',num2str(diagnostics.DKL.X,'%4.3e'),'\n ']);
end
if ~isnan(diagnostics.efficiency.X0)
    str{10} = sprintf([str{10},'    - initial conditions: ',num2str(diagnostics.efficiency.X0,'%4.3e'),'\n ']);
%     str{11} = sprintf([str{11},'    - initial conditions: ',num2str(diagnostics.DKL.X0,'%4.3e'),'\n ']);
end
if ~isnan(diagnostics.efficiency.Theta)
    str{10} = sprintf([str{10},'    - evolution parameters: ',num2str(diagnostics.efficiency.Theta,'%4.3e'),'\n ']);
%     str{11} = sprintf([str{11},'    - evolution parameters: ',num2str(diagnostics.DKL.Theta,'%4.3e'),'\n ']);
end
if ~isnan(diagnostics.efficiency.Phi)
    str{10} = sprintf([str{10},'    - observation parameters: ',num2str(diagnostics.efficiency.Phi,'%4.3e'),'\n ']);
%     str{11} = sprintf([str{11},'    - observation parameters: ', num2str(diagnostics.DKL.Phi,'%4.3e'),'\n ']);
end
if ~isnan(diagnostics.efficiency.alpha)
    str{10} = sprintf([str{10},'    - state noise precision hyperparameter: ',num2str(diagnostics.efficiency.alpha,'%4.3e'),'\n ']);
%     str{11} = sprintf([str{11},'    - state noise precision hyperparameter: ',num2str(diagnostics.DKL.alpha,'%4.3e'),'\n ']);
end
if ~isnan(diagnostics.efficiency.sigma)
    str{10} = sprintf([str{10},'    - data noise precision hyperparameter: ',num2str(diagnostics.efficiency.sigma','%4.3e '),'\n ']);
%     str{11} = sprintf([str{11},'    - data noise precision hyperparameter: ',num2str(diagnostics.DKL.sigma','%4.3e '),'\n ']);
end
str{11} = sprintf(['Classical fit accuracy metrics:','\n ']);
str{11} = sprintf([str{11},' - coefficient of determination (R2): ',num2str(out.fit.R2,'%4.3f'),'\n ']);
str{11} = sprintf([str{11},' - log-likelihood: ',num2str(out.fit.LL,'%4.3e'),'\n ']);
str{11} = sprintf([str{11},' - AIC: ',num2str(out.fit.AIC,'%4.3e'),'\n ']);
str{11} = sprintf([str{11},' - BIC: ',num2str(out.fit.BIC,'%4.3e'),'\n ']);
uicontrol('parent',hf,'style','text','tag','VBLaplace','units','normalized','position',[0.1,0.05,0.8,0.85],'backgroundcolor',[1,1,1],'HorizontalAlignment','left','fontsize',11,'string',str);


function myDeterministic(hfig)

try, hfig; catch, hfig = get(gco,'parent'); end
% first: clear diagnostics display
hc = intersect(findobj('tag','VBLaplace'),get(hfig,'children'));
if ~isempty(hc)
    delete(hc)
end

% Second: re-display VB-Laplace inversion output
ud = get(hfig,'userdata');
y = ud.out.y;
posterior = ud.out.options.init.posterior;
out = ud.out.options.init.out;

options = out.options;
options.noPause = 1;
options.DisplayWin = 1;
dim = out.dim;
suffStat = out.suffStat;
suffStat.dx0 = out.options.priors.muX0 - posterior.muX0;
suffStat.dtheta = out.options.priors.muTheta - posterior.muTheta;
suffStat.dphi = out.options.priors.muPhi - posterior.muPhi;
posterior.a_alpha = Inf;
posterior.b_alpha = 0;



% Initialize display figure
options.display.hfp = hfig;
options.figName = get(hfig,'name');
[options] = VBA_initDisplay(options);
delete(options.display.htt)
delete(options.display.hpause)
delete(options.display.hm)
delete(options.display.ho)
if options.dim.n == 0 || isinf(posterior.a_alpha(end))
    try delete(options.display.ha(8)); end
end
hfig = options.display.hfp;
drawnow

% Display data and hidden states (if any)
if options.dim.n > 0
    VBA_updateDisplay(posterior,suffStat,options,y,0,'X')
end
% Display precision hyperparameters
VBA_updateDisplay(posterior,suffStat,options,y,0,'precisions')
if ~options.OnLine && ~options.binomial
    xlabel(options.display.ha(6),' ')
    try
        xlabel(options.display.ha(8),' ')
    end
end
% Display model evidence
VBA_updateDisplay(posterior,suffStat,options,y,0,'F')
% Display parameters
if dim.n_theta >= 1
    VBA_updateDisplay(posterior,suffStat,options,y,0,'theta')
end
if dim.n_phi >= 1
    VBA_updateDisplay(posterior,suffStat,options,y,0,'phi')
end
try
    getSubplots
end

hdet = uicontrol('parent',hfig,'style','pushbutton','tag','VBLaplace','units','normalized','position',[0.4,0.05,0.2,0.02],'backgroundcolor',.8*[1,1,1],'string','diagnose deterministic?','callback',@diagnoseDeterministic);


function diagnoseDeterministic(ho,e)
ud = get(get(ho,'parent'),'userdata');
posterior = ud.out.options.init.posterior;
out = ud.out.options.init.out;
VBA_ReDisplay(posterior,out,1);


function myPriors(hfig)

try
    hfig;
catch
    hfig = get(gco,'parent');
end
% first: clear diagnostics display
hc = intersect(findobj('tag','VBLaplace'),get(hfig,'children'));
if ~isempty(hc)
    delete(hc)
end

% Second: re-display VB-Laplace inversion output
ud = get(hfig,'userdata');
out = ud.out;
y = out.y;
posterior = out.options.priors;
options = out.options;
options.noPause = 1;
options.DisplayWin =1;
dim = out.dim;
suffStat = out.suffStat;
suffStat.gx = ud.diagnostics.pgx;
% set dx = -prior.muX (for display purposes)
suffStat.dx0 = -posterior.muX0;
suffStat.dtheta = -posterior.muTheta;
suffStat.dphi = -posterior.muPhi;
suffStat.vy = ud.diagnostics.pvy;
try
    Ns=numel(options.sources);
catch
    Ns=1;
end

% Initialize display figure
options.display.hfp = hfig;
options.figName = get(hfig,'name');
[options] = VBA_initDisplay(options,1);
delete(options.display.htt)
delete(options.display.hpause)
delete(options.display.hm)
delete(options.display.ho)
if options.dim.n == 0 || isinf(posterior.a_alpha(end))
    try delete(options.display.ha(2*Ns+6)); end
end

% Display data and hidden states (if any)
if options.dim.n > 0
    options.OnLine = 0;
    VBA_updateDisplay(posterior,suffStat,options,y,0,'X')
end

% Display precision hyperparameters
VBA_updateDisplay(posterior,suffStat,options,y,0,'precisions')
if ~options.OnLine && ~options.binomial
    xlabel(options.display.ha(2*Ns+4),' ')
    try
        xlabel(options.display.ha(2*Ns+6),' ')
    end
end

% Display model evidence
VBA_updateDisplay(posterior,suffStat,options,y,0,'F')

% Display parameters
if dim.n_theta >= 1
    VBA_updateDisplay(posterior,suffStat,options,y,0,'theta')
end
if dim.n_phi >= 1
    VBA_updateDisplay(posterior,suffStat,options,y,0,'phi')
end

try
    getSubplots
end



function myConv(hfp)

try
    hf = hfp;
catch
    hf = get(gco,'parent');
end
hc = intersect(findobj('tag','VBLaplace'),get(hf,'children'));
if ~isempty(hc)
    delete(hc)
end

ud = get(hf,'userdata');
out = ud.out;
diagnostics = ud.diagnostics;

if length(out.suffStat.F)>2
    nit = length(out.suffStat.F)-1;
    ha = axes('parent',hf,'units','normalized','tag','VBLaplace','position',[0.15,0.6,0.5,0.3],'nextplot','add','xlim',[0,nit],'xtick',[0,nit],'xticklabel',{'prior (0)',['posterior (',num2str(nit),')']},'xgrid','off','ygrid','on');
    plot(ha,[0:nit],out.suffStat.F)
    plot(ha,[0:nit],out.suffStat.F,'.')
    [haf,hp1,hp2] = plotUncertainTimeSeries(diagnostics.LLH0*ones(1,2),3^2*ones(1,2),[0,nit],ha);
    set(hp1,'facecolor',[1 0 0])
    set(hp2,'color',[1 0 0])
    text(nit/2,diagnostics.LLH0-3/2,'log p(y|H0)','color',[1 0 0],'parent',ha);
    if ~out.options.OnLine
        title(ha,'VB optimization: F values','fontsize',11)
        xlabel(ha,'inner (Gauss-Newton) iterations')
    else
        title(ha,'online VB: F values','fontsize',11)
        xlabel(ha,'time samples')
    end
    ylabel(ha,'Free energy')
    box(ha,'off')
    if ~out.options.OnLine
        xtl = {'first iteration','last iteration'};
    else
        xtl = {'first time point','last time point'};
    end
    ha = axes('parent',hf,'units','normalized','tag','VBLaplace','position',[0.15,0.15,0.5,0.3],'nextplot','add','xlim',[0,nit-1],'xtick',[0,nit-1],'xticklabel',xtl,'xgrid','off','ygrid','on');
    plot(ha,[0:nit-1],diff(out.suffStat.F))
    plot(ha,[0:nit-1],diff(out.suffStat.F),'.')
    if ~out.options.OnLine
        title(ha,'VB optimization: F increments','fontsize',11)
        xlabel(ha,'inner (Gauss-Newton) iterations')
    else
        title(ha,'online VB: F increments','fontsize',11)
        xlabel(ha,'time samples')
    end
    ylabel(ha,'Free energy differences')
    box(ha,'off')
    try
        getSubplots
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
uicontrol('parent',hf,'style','text','tag','VBLaplace','units','normalized','position',[0.75,0.1,0.2,0.8],'backgroundcolor',[1,1,1],'HorizontalAlignment','left','fontsize',11,'string',str);



function myKernels()
% first: clear VB-Laplace inversion output display
hf = get(gco,'parent');
hc = intersect(findobj('tag','VBLaplace'),get(hf,'children'));
if ~isempty(hc)
    delete(hc)
end
ud = get(hf,'userdata');
dim = ud.out.options.dim;
unames = cell(dim.u,1);
for i=1:dim.u
    unames{i} = ['#',num2str(i)];
end
handles(1) = uicontrol('style','popupmenu','parent',hf,'tag','VBLaplace','units','normalized','position',[0.85 0.9 0.10 0.02],'fontsize',12,'string',unames,'callback',@myKerneli);
handles(2) = uicontrol('style','text','parent',hf,'tag','VBLaplace','BackgroundColor',get(hf,'color'),'units','normalized','position',[0.82 0.93 0.16 0.02],'fontsize',12,'string','display input...');
feval(@myKerneli,handles(1),[])



function myKerneli(hObject,evt)
hf = get(hObject,'parent');
ind = get(hObject,'Value');
ud = get(hf,'userdata');
try
    delete(ud.handles.hkernels)
end
out = ud.out;
kernels = ud.diagnostics.kernels;

% input effects - hidden states
kx = kernels.K1(:,:,ind);
handles.hkernels(1) = subplot(2,1,1,'parent',hf,'nextplot','add','ygrid','on','tag','VBLaplace');
plot(handles.hkernels(1),kernels.tgrid,kx)
set(handles.hkernels(1),'XLim',[min(kernels.tgrid) max(kernels.tgrid)])
pos = get(handles.hkernels(1),'position');
set(handles.hkernels(1),'position',[0.2 pos(2) 0.6 pos(4)])
title(handles.hkernels(1),['states impulse responses to input #''' num2str(ind) ''''],'fontsize',12)
xlabel(handles.hkernels(1),'time')

% input effects - observables
ky = kernels.H1(:,:,ind);
handles.hkernels(2) = subplot(2,1,2,'parent',hf,'nextplot','add','ygrid','on','tag','VBLaplace');
pos = get(handles.hkernels(2),'position');
set(handles.hkernels(2),'position',[0.2 pos(2) 0.6 pos(4)])
plot(handles.hkernels(2),kernels.tgrid,ky)
set(handles.hkernels(2),'XLim',[min(kernels.tgrid) max(kernels.tgrid)])
title(handles.hkernels(2),['observables impulse responses to input #''' num2str(ind) ''''],'fontsize',12)
xlabel(handles.hkernels(2),'time')

ud.handles = handles;
set(hf,'userdata',ud);
try getSubplots; end



function myDiagnostics()

% first: clear VB-Laplace inversion output display
hf = get(gco,'parent');
hc = intersect(findobj('tag','VBLaplace'),get(hf,'children'));
if ~isempty(hc)
    delete(hc)
end

% Second: display diagnostics
ud = get(hf,'userdata');
out = ud.out;
y = out.y;
diagnostics = ud.diagnostics;

% display micro-time hidden-states
if ~isempty(diagnostics.MT_x)
    display.ha(1) = subplot(4,2,1,'parent',hf,'nextplot','add','tag','VBLaplace','ygrid','on','box','off');
    title(display.ha(1),'micro-time resolution predicted data','fontsize',11)
    xlabel(display.ha(1),'time','fontsize',8)
    ylabel(display.ha(1),'g(x) & y','fontsize',8)
    plot(display.ha(1),diagnostics.microTime,diagnostics.MT_gx')
    plot(display.ha(1),diagnostics.microTime(diagnostics.sampleInd),diagnostics.MT_gx(:,diagnostics.sampleInd)','.')
    plot(display.ha(1),diagnostics.microTime(diagnostics.sampleInd),y,':')
    axis(display.ha(1),'tight')
    display.ha(2) = subplot(4,2,2,'parent',hf,'nextplot','add','tag','VBLaplace','ygrid','on','box','off');
    title(display.ha(2),'micro-time resolution hidden states','fontsize',11)
    xlabel(display.ha(2),'time','fontsize',8)
    ylabel(display.ha(2),'x','fontsize',8)
    plot(display.ha(2),diagnostics.microTime,diagnostics.MT_x')
    plot(display.ha(2),diagnostics.microTime(diagnostics.sampleInd),diagnostics.MT_x(:,diagnostics.sampleInd)','.')
    axis(display.ha(2),'tight')
end

% display data noise
xlim = [diagnostics.dy.nx(1)-diagnostics.dy.d,diagnostics.dy.nx(end)+diagnostics.dy.d];
display.ha(3) = subplot(4,2,5,'parent',hf,'nextplot','add','xlim',xlim,'ygrid','on','tag','VBLaplace','box','off');
title(display.ha(3),'residuals empirical distribution','fontsize',11)
xlabel(display.ha(3),'e(t) = y(t)-g(x(t))','fontsize',8)
ylabel(display.ha(3),'p(e|y)','fontsize',8)
bar(diagnostics.dy.nx,diagnostics.dy.ny,'facecolor',[.8 .8 .8],'parent',display.ha(3))
plot(display.ha(3),diagnostics.dy.grid,diagnostics.dy.pg,'r')
if ~out.options.binomial
    plot(display.ha(3),diagnostics.dy.grid,diagnostics.dy.pg2,'g')
end
if ~out.options.binomial
    legend(display.ha(3),{'empirical histogram','Gaussian approx','posterior approx'})
else
    legend(display.ha(3),{'empirical histogram','Gaussian approx'})
end

if out.options.dim.n > 0
    gri = diagnostics.microTime(diagnostics.sampleInd);
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
display.ha(7) = subplot(4,2,3,'parent',hf,'nextplot','add','tag','VBLaplace','ygrid','on','box','off');
plot(display.ha(7),gri,out.suffStat.dy')
axis(display.ha(7),'tight')
title(display.ha(7),'residuals time series','fontsize',11)
xlabel(display.ha(7),ti,'fontsize',8)
ylabel(display.ha(7),'e(t) = y(t)-g(x(t))','fontsize',8)

% display autocorrelation of residuals
if ~isweird(diagnostics.dy.R) && out.dim.n_t > 1
    display.ha(5) = subplot(4,2,7,'parent',hf);
    plot(display.ha(5),[-out.options.dim.n_t:out.options.dim.n_t-1],fftshift(diagnostics.dy.R)')
    axis(display.ha(5),'tight')
    title(display.ha(5),'residuals empirical autocorrelation','fontsize',11)
    xlabel(display.ha(5),'lag tau','fontsize',8)
    ylabel(display.ha(5),'Corr[e(t),e(t+tau)]','fontsize',8)
    set(display.ha(5),...
        'tag','VBLaplace',...
        'ygrid','on',...
        'box','off');
end


% display state noise
if ~isempty(diagnostics.dx.dx)
    xlim = [diagnostics.dx.nx(1)-diagnostics.dx.d,diagnostics.dx.nx(end)+diagnostics.dx.d];
    display.ha(4) = subplot(4,2,6,'parent',hf,'nextplot','add','xlim',xlim,'ygrid','on','tag','VBLaplace','box','off');
    title(display.ha(4),'state noise empirical distribution','fontsize',11)
    xlabel(display.ha(4),'eta(t) = x(t+1)-f(x(t))','fontsize',8)
    ylabel(display.ha(4),'p(eta|y)','fontsize',8)
    bar(diagnostics.dx.nx,diagnostics.dx.ny,'facecolor',[.8 .8 .8],'parent',display.ha(4))
    plot(display.ha(4),diagnostics.dx.grid,diagnostics.dx.pg,'r')
    plot(display.ha(4),diagnostics.dx.grid,diagnostics.dx.pg2,'g')
    legend(display.ha(4),{'empirical histogram','Gaussian approx','posterior approx'})
    
    display.ha(8) = subplot(4,2,4,'parent',hf,'nextplot','add','tag','VBLaplace','ygrid','on','box','off');
    try
        plotUncertainTimeSeries(out.suffStat.dx,out.suffStat.vdx,diagnostics.microTime(diagnostics.sampleInd),display.ha(8));
        
    catch
        plot(display.ha(8),diagnostics.microTime(diagnostics.sampleInd),out.suffStat.dx')
    end
    axis(display.ha(8),'tight')
    title(display.ha(8),'state noise time series','fontsize',11)
    xlabel(display.ha(8),'time','fontsize',8)
    ylabel(display.ha(8),'eta(t) = x(t+1)-f(x(t))','fontsize',8)
    
end


% display parameters posterior correlation matrix
display.ha(6) = subplot(4,2,8,'parent',hf);
imagesc(diagnostics.C,'parent',display.ha(6))
title(display.ha(6),'parameters posterior correlation matrix','fontsize',11)
set(display.ha(6),...
    'tag','VBLaplace',...
    'xtick',diagnostics.ltick,...
    'ytick',diagnostics.ltick,...
    'xticklabel',diagnostics.ticklabel,...
    'yticklabel',diagnostics.ticklabel,...
    'box','off',...
    'nextplot','add');
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

try
    getSubplots
end


function myVB(hfig)

try
    hfig;
catch
    hfig = get(gco,'parent');
end
% first: clear diagnostics display
hc = intersect(findobj('tag','VBLaplace'),get(hfig,'children'));
if ~isempty(hc)
    delete(hc)
end


% Second: re-display VB-Laplace inversion output
ud = get(hfig,'userdata');
out = ud.out;
y = out.y;
posterior = ud.posterior;

options = out.options;
options.noPause = 1;
options.DisplayWin =1;
suffStat = out.suffStat;
try
    Ns=numel(options.sources);
catch
    Ns=1;
end
dim = out.dim;

% Initialize display figure
options.display.hfp = hfig;
options.figName = get(hfig,'name');
[options] = VBA_initDisplay(options);
delete(options.display.htt)
delete(options.display.hpause)
delete(options.display.hm)
delete(options.display.ho)
if options.dim.n == 0 || isinf(posterior.a_alpha(end))
    try delete(options.display.ha(2*Ns+6)); end
end
hfig = options.display.hfp;
drawnow

% Display data and hidden states (if any)
if options.dim.n > 0
    VBA_updateDisplay(posterior,suffStat,options,y,0,'X')
end

% Display precision hyperparameters
VBA_updateDisplay(posterior,suffStat,options,y,0,'precisions')
if ~options.OnLine && ~options.binomial
    xlabel(options.display.ha(2*Ns+4),' ')
    try
        xlabel(options.display.ha(2*Ns+6),' ')
    end
end

% Display model evidence
VBA_updateDisplay(posterior,suffStat,options,y,0,'F')

% Display parameters
if dim.n_theta >= 1
    VBA_updateDisplay(posterior,suffStat,options,y,0,'theta')
end
if dim.n_phi >= 1
    VBA_updateDisplay(posterior,suffStat,options,y,0,'phi')
end

try
    getSubplots
end



function diagnostics = getDiagnostics(posterior,out)

if out.options.verbose
    fprintf(1,['Deriving diagnostics ...'])
end

u = out.u;
y = out.y;

try; out.fit; catch; out.fit = VBA_fit(posterior,out); end
    
% get kernels (NB: dcm = special case)
if isequal(out.options.f_fname,@f_DCMwHRF) && isequal(out.options.g_fname,@g_HRF3)
    dcm = 1;
    [out.options] = VBA_check4DCM(out.options);
else
    dcm = 0;
end
[kernels.H1,kernels.K1,kernels.tgrid] = getKernels(posterior,out,dcm);

% get null model (H0) evidence
[LLH0] = VBA_LMEH0(y,out.options);

% Entropies and KL divergences

if sum([out.options.sources(:).type]==0)>0 %~out.options.binomial
    efficiency.sigma = -out.suffStat.Ssigma;
    m0 = out.options.priors.a_sigma./out.options.priors.b_sigma;
    v0 = out.options.priors.a_sigma./out.options.priors.b_sigma.^2;
    m = posterior.a_sigma./posterior.b_sigma;
    v = posterior.a_sigma./posterior.b_sigma.^2;
    DKL.sigma = VB_KL(m0,v0,m,v,'Gamma');
else
    efficiency.sigma = NaN;
    DKL.sigma = NaN;
end

if out.dim.n > 0 % hidden states and initial conditions
    efficiency.X = -out.suffStat.SX;
    efficiency.X0 = -out.suffStat.SX0;
    if isinf(out.options.priors.a_alpha) ...
            && isequal(out.options.priors.b_alpha,0)
        efficiency.alpha = NaN;
        DKL.alpha = NaN;
    else
        efficiency.alpha = -out.suffStat.Salpha;
        m0 = out.options.priors.a_alpha./out.options.priors.b_alpha;
        v0 = out.options.priors.a_alpha./out.options.priors.b_alpha^2;
        m = posterior.a_alpha(end)./posterior.b_alpha(end);
        v = posterior.a_alpha(end)./posterior.b_alpha(end)^2;
        DKL.alpha = VB_KL(m0,v0,m,v,'Gamma');
    end
    try
        DKL.X = 0;
        for t=1:out.dim.n_t
            IN = out.options.params2update.x{t};
            m0 = out.options.priors.muX(IN,t);
            v0 = out.options.priors.SigmaX.current{t}(IN,IN);
            m = posterior.muX(IN,t);
            v = posterior.SigmaX.current{t}(IN,IN);
            DKL.X = DKL.X + VB_KL(m0,v0,m,v,'Normal');
        end
    catch
        DKL.X = NaN;
    end
    IN = out.options.params2update.x0;
    m0 = out.options.priors.muX0(IN);
    v0 = out.options.priors.SigmaX0(IN,IN);
    m = posterior.muX0(IN);
    v = posterior.SigmaX0(IN,IN);
    DKL.X0 = VB_KL(m0,v0,m,v,'Normal');
else
    efficiency.X = NaN;
    efficiency.X0 = NaN;
    efficiency.alpha = NaN;
    DKL.X = NaN;
    DKL.X0 = NaN;
    DKL.alpha = NaN;
end
if out.dim.n_phi > 0 % observation parameters
    efficiency.Phi = -out.suffStat.Sphi;
    IN = out.options.params2update.phi;
    m0 = out.options.priors.muPhi(IN);
    v0 = out.options.priors.SigmaPhi(IN,IN);
    if ~out.options.OnLine
        m = posterior.muPhi(IN);
        v = posterior.SigmaPhi(IN,IN);
    else
        m = posterior.muPhi(IN,end);
        v = posterior.SigmaPhi{end}(IN,IN);
    end
    DKL.Phi = VB_KL(m0,v0,m,v,'Normal');
else
    efficiency.Phi = NaN;
    DKL.Phi = NaN;
end
if out.dim.n_theta > 0 % evolution parameters
    efficiency.Theta = -out.suffStat.Stheta;
    IN = out.options.params2update.theta;
    m0 = out.options.priors.muTheta(IN);
    v0 = out.options.priors.SigmaTheta(IN,IN);
    if ~out.options.OnLine
        m = posterior.muTheta(IN);
        v = posterior.SigmaTheta(IN,IN);
    else
        m = posterior.muTheta(IN,end);
        v = posterior.SigmaTheta{end}(IN,IN);
    end
    DKL.Theta = VB_KL(m0,v0,m,v,'Normal');
else
    efficiency.Theta = NaN;
    DKL.Theta = NaN;
end

% get prior predictive density
[muy,Vy] = VBA_getLaplace(u,out.options.f_fname,out.options.g_fname,out.dim,out.options);

% get micro-time posterior hidden-states estimates
try
    [MT_x,MT_gx,microTime,sampleInd] = VBA_microTime(posterior,u,out);
catch
    MT_x = [];
    MT_gx = [];
    microTime = [];
    sampleInd = [];
end

% get residuals: data noise
dy.dy = out.suffStat.dy(:);
dy.R = spm_autocorr(out.suffStat.dy);
dy.m = mean(dy.dy);
dy.v = var(dy.dy);
[dy.ny,dy.nx] = hist(dy.dy,10);
dy.ny = dy.ny./sum(dy.ny);
d = diff(dy.nx);
d = abs(d(1));
dy.d = d;
spgy = sum(exp(-0.5.*(dy.m-dy.nx).^2./dy.v));
dy.grid = dy.nx(1):d*1e-2:dy.nx(end);
dy.pg = exp(-0.5.*(dy.m-dy.grid).^2./dy.v);
dy.pg = dy.pg./spgy;
if  ~out.options.binomial
    shat = posterior.a_sigma(end)./posterior.b_sigma(end);
    spgy = sum(exp(-0.5.*shat.*dy.nx.^2));
    dy.pg2 = exp(-0.5.*shat.*dy.grid.^2);
    dy.pg2 = dy.pg2./spgy;
end

% get residuals: state noise
dx.dx = out.suffStat.dx(:);
if ~isempty(dx.dx)
    dx.m = mean(dx.dx);
    dx.v = var(dx.dx);
    [dx.ny,dx.nx] = hist(dx.dx,10);
    dx.ny = dx.ny./sum(dx.ny);
    d = diff(dx.nx);
    d = abs(d(1));
    dx.d = d;
    spgy = sum(exp(-0.5.*(dx.m-dx.nx).^2./dx.v));
    dx.grid = dx.nx(1):d*1e-2:dx.nx(end);
    dx.pg = exp(-0.5.*(dx.m-dx.grid).^2./dx.v);
    dx.pg = dx.pg./spgy;
    ahat = posterior.a_alpha(end)./posterior.b_alpha(end);
    spgy = sum(exp(-0.5.*ahat.*dx.nx.^2));
    dx.pg2 = exp(-0.5.*ahat.*dx.grid.^2);
    dx.pg2 = dx.pg2./spgy;
end


% get parameters posterior correlation matrix
if out.dim.n > 0 && isinf(out.options.priors.a_alpha) && isequal(out.options.priors.b_alpha,0)
    S = out.suffStat.ODE_posterior.SigmaPhi;
else
    S = NaN*zeros(out.dim.n+out.dim.n_theta+out.dim.n_phi);
    ind = 0;
    if out.dim.n_phi > 0
        if iscell(posterior.SigmaPhi) % online version
            SP = posterior.SigmaPhi{end};
        else
            SP = posterior.SigmaPhi;
        end
        S(1:out.dim.n_phi,1:out.dim.n_phi) = SP;
        ind = out.dim.n_phi;
    end
    if out.dim.n_theta > 0
        if iscell(posterior.SigmaTheta) % online version
            SP = posterior.SigmaTheta{end};
        else
            SP = posterior.SigmaTheta;
        end
        S(ind+1:ind+out.dim.n_theta,ind+1:ind+out.dim.n_theta) = SP;
        ind = ind + out.dim.n_theta;
    end
    if out.dim.n > 0 && out.options.updateX0
        if iscell(posterior.SigmaX0) % online version
            SP = posterior.SigmaX0{end};
        else
            SP = posterior.SigmaX0;
        end
        S(ind+1:ind+out.dim.n,ind+1:ind+out.dim.n) = SP;
    end
end
C = cov2corr(S);
C = C + diag(NaN.*diag(C));
tick = [0];
ltick = [];
ticklabel = cell(0,0);
if out.dim.n_phi > 0
    ltick = [ltick,tick(end)+out.dim.n_phi/2];
    tick = [tick,out.dim.n_phi];
    ticklabel{end+1} = 'phi';
end
if out.dim.n_theta > 0
    ltick = [ltick,tick(end)+out.dim.n_theta/2];
    tick = [tick,tick(end)+out.dim.n_theta];
    ticklabel{end+1} = 'theta';
end
if out.dim.n > 0 && out.options.updateX0
    ltick = [ltick,tick(end)+out.dim.n/2];
    tick = [tick,tick(end)+out.dim.n];
    ticklabel{end+1} = 'x0';
end
tick = tick +0.5;
tick = tick(2:end-1);
ltick = ltick + 0.5;

diagnostics.pgx = reshape(muy,out.dim.p,out.dim.n_t);
diagnostics.pvy = reshape(diag(Vy),out.dim.p,out.dim.n_t);
if ~isempty(kernels.tgrid)
    diagnostics.kernels = kernels;
else
    diagnostics.kernels = [];
end
diagnostics.efficiency = efficiency;
diagnostics.DKL = DKL;
diagnostics.LLH0 = LLH0;
diagnostics.MT_x = MT_x;
diagnostics.MT_gx = MT_gx;
diagnostics.microTime = microTime;
diagnostics.sampleInd = sampleInd;
diagnostics.dy = dy;
diagnostics.dx = dx;
diagnostics.ltick = ltick;
diagnostics.tick = tick;
diagnostics.ticklabel = ticklabel;
diagnostics.C = C;

if out.options.verbose
    fprintf(' OK.')
    fprintf('\n')
end


