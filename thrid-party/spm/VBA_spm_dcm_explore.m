function [] = VBA_spm_dcm_explore(P)

% build GUI for DCM results reviewing

if ~nargin
    [P, sts] = VBA_spm_select(Inf,'^DCM.*\.mat$','select DCM_???.mat');
    if ~sts, return; end
    nf = size(P,1);
    P0 = P;
    P = cell(nf,1);
    for i=1:nf
        [pathstr,name,ext] = fileparts(P0(i,:));
        P{i} = [pathstr,filesep,name,'.mat'];
    end
end
if ~isstruct(P)
    if iscell(P)
        nf = length(P);
        for i=1:nf
            VBA_spm_dcm_explore(P{i})
        end
        return
    else
        try
            load(P);
            DCM.filename = P;
        catch
            fprintf('Error: Cannot load DCM file "%s".',P)
            disp(' ')
            return
        end
    end
else
    DCM = P;
%     DCM.filename = [];
end


% Get model specification structure (see spm_nlsi)
%--------------------------------------------------------------------------

try
    ud.dim.m  = DCM.M.m; % number of inputs
    ud.dim.n  = DCM.M.n; % number of hidden states
    ud.dim.l  = DCM.M.l; % number of regions (responses)
    ud.estimated = 1;
catch
    DCM = fillinDcm(DCM);
    ud.dim.m  = DCM.M.m; % number of inputs
    ud.dim.n  = DCM.M.n; % number of hidden states
    ud.dim.l  = DCM.M.l; % number of regions (responses)
    ud.estimated = 0;
end
ud.display.U  = 0.9; % p-value threshold for display
if ~isfield(DCM,'filename')
    DCM.filename = '[ DCM worskpace variable ]';
end
if ~isfield(DCM,'d')
    DCM.d = [];
end
try
    if DCM.options.two_state
        ud.Ep.A = exp(DCM.Ep.A);
        ud.Ep.B = exp(DCM.Ep.B);
        ud.Ep.D = exp(DCM.Ep.D);
        ud.Ep.C = DCM.Ep.C;
        disp('NB: The (A,B,D) parameters are scale parameters')
    else
        ud.Ep.A = DCM.Ep.A;
        ud.Ep.B = DCM.Ep.B;
        ud.Ep.D = DCM.Ep.D;
        ud.Ep.C = DCM.Ep.C;
    end
catch % old DCM structure
    ud.Ep.A = DCM.A;
    ud.Ep.B = DCM.B;
    ud.Ep.C = DCM.C;
    DCM.Pp.A = DCM.pA;
    DCM.Pp.B = DCM.pB;
    DCM.Pp.C = DCM.pC;
    try
        ud.Ep.D = DCM.D;
        DCM.Pp.D = DCM.pD;
    catch
        ud.Ep.D = [];
        DCM.Pp.D = [];
        DCM.d = zeros(ud.dim.l,ud.dim.l,ud.dim.l);
    end
    DCM.options.stochastic = 0;
    DCM.options.two_state = 0;
    try
        DCM.TE;
    catch
        DCM.TE = 0;
    end
    try
        DCM.M.IS;
    catch
        DCM.M.IS = '?';
    end
end

% test DCM parameters
try
    DCM = VBA_spm_dcm_test(DCM);
catch
    DCM = myspm_dcm_test(DCM);
end
ud.Pp = DCM.sdr; % posterior proba for P=0 (using SD-ratios)
ud.Pcorr = VBA_cov2corr(DCM.Cp); % posterior correlation matrix
for i=1:ud.dim.l
    PSS(i)  = sum(DCM.y(:,i).^2);
    RSS(i)  = sum(DCM.R(:,i).^2);
    ud.R2.reg(i) = PSS(i)/(PSS(i) + RSS(i)); % coef of determination
end
ud.R2.tot =  sum(PSS)/sum(PSS+RSS);

pos0 = get(0,'screenSize');
pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];
[pathstr,filename,ext] = fileparts(DCM.filename);
figname = ['Explore ',filename];
hfp = figure('position',pos,'color',[1 1 1],'name',figname,'tag','exploreDCM','Renderer','OpenGL','menu','none');
uimenu(hfp,'Label','Open new DCM file','Callback','VBA_spm_dcm_explore');
if isfield(DCM,'VBA')
    uimenu(hfp,'Label','Diagnose VBA','Callback',@vba_display);
end
if ud.estimated
    labels = {'summary','parameters','kernels','diagnostics'};
    callbacks = {@mySummary,@myNeuro,@myHemo,@myDiagnostics};
else
    labels = {'summary','parameters'};
    callbacks = {@mySummary,@myNeuro};
end
[ud.tabs.handles] = VBA_spm_uitab(hfp,labels,callbacks,'exploreDCM',1);
set(ud.tabs.handles.htab,'backgroundcolor',[1 1 1])
set(ud.tabs.handles.hh,'backgroundcolor',[1 1 1])
set(ud.tabs.handles.hp,'HighlightColor',0.8*[1 1 1])
set(ud.tabs.handles.hp,'backgroundcolor',[1 1 1])
ud.DCM = DCM;
set(hfp,'userdata',ud);
feval(@mySummary,hfp)


function vba_display(hObject,evt)
ud = get(get(hObject,'parent'),'userdata');
VBA_ReDisplay(ud.DCM.VBA.posterior,ud.DCM.VBA.out,1)


function mySummary(hfp)
try
    hf = hfp;
catch
    hf = get(gco,'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles.hp,'children'))
labels = {'spec','ROIs','inputs','outputs'};
callbacks = {@mySpec,@myROIs,@myI,@myO};
[ud.tabs.handles2] = VBA_spm_uitab(ud.tabs.handles.hp,labels,callbacks,'exploreDCM0',1);
set(ud.tabs.handles2.htab,'backgroundcolor',[1 1 1])
set(ud.tabs.handles2.hh,'backgroundcolor',[1 1 1])
set(ud.tabs.handles2.hp,'HighlightColor',0.8*[1 1 1])
set(ud.tabs.handles2.hp,'backgroundcolor',[1 1 1])
set(hf,'userdata',ud);
feval(@mySpec,hf)



function myNeuro(hfp)

try
    hf = hfp;
catch
    hf = get(gco,'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles.hp,'children'))
if ud.estimated
    labels = {'A','B','C','D'};
    callbacks = {@myA,@myB,@myC,@myD};
else
    labels = {'A','B','C','D'};
    callbacks = {@myA,@myB,@myC,@myD};
end
[ud.tabs.handles2] = VBA_spm_uitab(ud.tabs.handles.hp,labels,callbacks,'exploreDCM0',1);
set(ud.tabs.handles2.htab,'backgroundcolor',[1 1 1])
set(ud.tabs.handles2.hh,'backgroundcolor',[1 1 1])
set(ud.tabs.handles2.hp,'HighlightColor',0.8*[1 1 1])
set(ud.tabs.handles2.hp,'backgroundcolor',[1 1 1])
set(hf,'userdata',ud);
feval(@myA,hf)


function myHemo(hfp)
try
    hf = hfp;
catch
    hf = get(gco,'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles.hp,'children'))
DCM = ud.DCM;
nu= size(DCM.U.u,1);
if nu >= 1
    ud.tabs.handles3(1) = uicontrol(...
        'style','popupmenu',...
        'parent',ud.tabs.handles.hp,...
        'units','normalized',...
        'position',[0.85 0.9 0.12 0.02],...
        'fontsize',12,...
        'string',DCM.U.name,...
        'callback',@myKernel);
    ud.tabs.handles3(2) = uicontrol(...
        'style','text',...
        'parent',ud.tabs.handles.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.82 0.93 0.18 0.02],...
        'fontsize',12,...
        'string','display input...');
    set(hf,'userdata',ud)
    feval(@myKernel,ud.tabs.handles3(1),[])
else
    uicontrol(...
        'style','text',...
        'parent',ud.tabs.handles.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.25 0.5 0.5 0.1],...
        'fontsize',12,...
        'string','no input, no Volterra kernel!');
end


function myDiagnostics(hfp)
try
    hf = hfp;
catch
    hf = get(gco,'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles.hp,'children'))
labels = {'cov','residuals'};
callbacks = {@myCov,@myResiduals};
[ud.tabs.handles2] = VBA_spm_uitab(ud.tabs.handles.hp,labels,callbacks,'exploreDCM0',1);
set(ud.tabs.handles2.htab,'backgroundcolor',[1 1 1])
set(ud.tabs.handles2.hh,'backgroundcolor',[1 1 1])
set(ud.tabs.handles2.hp,'HighlightColor',0.8*[1 1 1])
set(ud.tabs.handles2.hp,'backgroundcolor',[1 1 1])
set(hf,'userdata',ud);
feval(@myCov,hf)



function myKernel(hObject,evt)
hf = get(get(hObject,'parent'),'parent');
ind = get(hObject,'Value');
ud = get(hf,'userdata');
delete(setdiff(get(ud.tabs.handles.hp,'children'),ud.tabs.handles3))
DCM = ud.DCM;

x     = [1:DCM.M.N] * DCM.M.dt;
d     = 2 / DCM.M.dt;

% input effects - neuronal
%--------------------------------------------------------------
y = DCM.K1(:,:,ind);
ha = subplot(2,1,1,'parent',ud.tabs.handles.hp);
plot(ha,x,y)
set(ha,'XLim',[1 max(x)])
pos = get(ha,'position');
set(ha,'position',[0.2 pos(2) 0.6 pos(4)])
title(ha,['neuronal impulse responses to input ''' DCM.U.name{ind} ''''],'fontsize',12)
grid(ha,'on')
xlabel(ha,'time (seconds)')
legend(ha,DCM.Y.name)

% input effects - hemodynamic
%--------------------------------------------------------------
y = DCM.K1(:,:,ind);
k = DCM.H1(:,:,ind);
ha = subplot(2,1,2,'parent',ud.tabs.handles.hp,'nextplot','add');
pos = get(ha,'position');
set(ha,'position',[0.2 pos(2) 0.6 pos(4)])
plot(ha,x,k)
plot(ha,x,y,':')
set(ha,'XLim',[1 max(x)])
title(ha,['hemodynamic impulse responses to input ''' DCM.U.name{ind} ''''],'fontsize',12)
grid(ha,'on')
xlabel(ha,'time (seconds)')
legend(ha,DCM.Y.name)
try;VBA_getSubplots; end



function mySpec(hfp)
try
    hf = hfp;
catch
    hf = get(get(gco,'parent'),'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles2.hp,'children'))
DCM = ud.DCM;
dc = diag(DCM.Cp);
np = length(find(dc~=0));

try
    str{1} = sprintf(['Date: ',datestr(DCM.date),'\n ']);
catch
    str{1} = '';
end
try
    [pathstr,filename,ext] = fileparts(DCM.filename);
    str{2} = sprintf(['File name: ',filename,ext,'\n ']);
catch
    str{2} = '';
end
str{3} = sprintf([...
    'Dimensions of the model:','\n ',...
    '    - data: p=',num2str(DCM.M.l),'\n ',...
    '    - time samples: t=',num2str(size(DCM.y,1)),'\n ',...
    '    - hidden states: n=',num2str(DCM.M.n),'\n ',...
    '    - inputs: nu=',num2str(DCM.M.m),'\n ',...
    '    - parameters: n_theta=',num2str(np),'\n ']);
str{4} = sprintf(['EPI repetition time: TR=',num2str(DCM.Y.dt),'\n',...
    'EPI echo time: TE=',num2str(DCM.TE),'\n']);
if ~DCM.options.stochastic
    str{5} = sprintf('This is a deterministic DCM');
else
    str{5} = sprintf('This is a stochastic DCM');
end
if DCM.options.two_state
    str{5} = sprintf([str{5},' (two states per region):']);
else
    str{5} = sprintf([str{5},' (one state per region):']);
end

if ud.estimated
    gfn = DCM.M.g;
    ffn = DCM.M.f;
    ifn = DCM.M.IS;
    str{6} = sprintf([...
        '    - observation function: ',gfn,'\n',...
        '    - evolution function: ',ffn,'\n',...
        '    - integration function: ',ifn,'\n ']);
else
    str{6} = [];
end
na = sum(DCM.a(:)) - sum(diag(DCM.a));
nb = sum(DCM.b(:));
nc = sum(DCM.c(:));
nd = sum(DCM.d(:));
str{7} = sprintf([...
    'DCM structure:','\n ',...
    '    - number of connections: ',num2str(na),' (+',num2str(DCM.M.l),' inhibitory self-connections)','\n ',...
    '    - number of modulatory effects: ',num2str(nb),'\n ',...
    '    - number of input-state couplings: ',num2str(nc),'\n ',...
    '    - number of (nonlinear) gating effects: ',num2str(nd),'\n ']);
str{8} = ['Variance explained:  R2=',num2str(100*ud.R2.tot,'%4.1f'),'%'];
if isnumeric(DCM.F)
    str{9} = ['Log model evidence: log p(y|m) > ',num2str(DCM.F,'%4.3e')];
else
    try
        str{9} = ['Warning: ',DCM.F];
    catch
        str{9} = ' ';
    end
end
uicontrol(...
    'parent',ud.tabs.handles2.hp,...
    'style','text',...
    'units','normalized',...
    'position',[0.1,0.25,0.8,0.68],...
    'backgroundcolor',[1,1,1],...
    'HorizontalAlignment','left',...
    'fontsize',12,...
    'string',str);


function myROIs(hfp)
try
    hf = hfp;
catch
    hf = get(get(gco,'parent'),'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles2.hp,'children'))
DCM = ud.DCM;
if isfield(DCM,'xY')
    ha(1) = subplot(2,1,1,...
        'parent',ud.tabs.handles2.hp,...
        'visible','off');
    set(ha,'position',get(ha,'position')+[0 -0.1 0 0.1])
    h = myspm_dcm_graph(DCM.xY,ud.Ep.A,ha);
    try; set(hf,'renderer','opengl'); end
    try
        set(h.handles.BUTTONS.transp,...
            'parent',ud.tabs.handles2.hp,...
            'position',get(h.handles.BUTTONS.transp,'position')-0.04*[1 1 0 0])
    end
    colormap(ha(1),'gray')
    hc = get(ha(1),'children');
    if isempty(hc)
        delete(ha)
    else
        set(ha,'visible','on')
        xlabel(ha,'x (mm)')
        ylabel(ha,'y (mm)')
        zlabel(ha,'z (mm)')
    end
    ha(2) = subplot(2,1,2,...
        'parent',ud.tabs.handles2.hp);
    set(ha(2),'position',get(ha(2),'position')+[0 -0.075 0 0])
    line([0 4],[0 0],'parent',ha(2))
    text(0,-1,'Name',...
        'FontSize',14,...
        'parent',ha(2))
    text(2,-1,'Voxels',...
        'FontSize',14,...
        'parent',ha(2))
    text(3,-1,'Location (mm)',...
        'FontSize',14,...
        'parent',ha(2))
    line([0 4],[-2 -2],...
        'LineWidth',4,...
        'parent',ha(2))
    y = -3;
    for i = 1:length(DCM.xY)
        name = DCM.xY(i).name;
        N    = length(DCM.xY(i).s);
        L    = DCM.xY(i).xyz;
        r    = DCM.xY(i).spec;
        text(0,y,name,...
            'FontWeight','bold',...
            'FontSize',12,...
            'parent',ha(2))
        text(2,y,sprintf('%0.0f',N),...
            'FontSize',12,...
            'parent',ha(2))
        text(3,y,sprintf('%-4.0f %-4.0f %-4.0f',L),...
            'FontSize',12,...
            'parent',ha(2))
        y = y - 1;
    end
    line([0 4],[y y],'parent',ha(2))
    axis(ha(2),'off')
try;VBA_getSubplots; end
else
    ht = uicontrol(...
        'style','text',...
        'parent',ud.tabs.handles2.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.25 0.5 0.5 0.1],...
        'fontsize',12,...
        'string','no information on ROI selection ');
end

function myI(hfp)
try
    hf = hfp;
catch
    hf = get(get(gco,'parent'),'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles2.hp,'children'))
DCM = ud.DCM;
%-priors-----------------------------------------------------------
x     = [1:length(DCM.U.u)]*DCM.U.dt;
t     = [1:length(DCM.Y.y)]*DCM.Y.dt;
for i = 1:size(DCM.U.u,2)
    ha = subplot(size(DCM.U.u,2),1,i,...
        'parent',ud.tabs.handles2.hp,...
        'nextplot','add');
    plot(ha,x,full(DCM.U.u(:,i)))
    if DCM.options.stochastic
        try
            plot(ha,t,DCM.qU.v{2}(i,:),':')
            legend(ha,{'prior','posterior'})
        end
    end
    title(ha,DCM.U.name{i},'FontSize',12)
    ylabel(ha,'event density (Hz)')
    xlabel(ha,'time (seconds)')
    a = axis(ha);
    axis(ha,[0 max(x) (a(3) - 1/8) (a(4) + 1/8)]);
    grid(ha,'on')
    box(ha,'on')
end
try;VBA_getSubplots; end



function myO(hfp)
try
    hf = hfp;
catch
    hf = get(get(gco,'parent'),'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles2.hp,'children'))
DCM = ud.DCM;
ud.tabs.handles3(1) = uicontrol('style','popupmenu',...
    'parent',ud.tabs.handles2.hp,...
    'units','normalized',...
    'fontsize',12,...
    'string',DCM.Y.name,...
    'callback',@myOi);
set(ud.tabs.handles3(1),'position',[0.85 0.9 0.12 0.02])
ud.tabs.handles3(2) = uicontrol('style','text',...
    'parent',ud.tabs.handles2.hp,...
    'BackgroundColor',get(hf,'color'),...
    'units','normalized',...
    'position',[0.82 0.93 0.18 0.02],...
    'fontsize',12,...
    'string','display ROI...');
set(hf,'userdata',ud)
feval(@myOi,ud.tabs.handles3(1),[])



function myOi(hObject,evt)
hf = get(get(get(hObject,'parent'),'parent'),'parent');
i = get(hObject,'Value');
ud = get(hf,'userdata');
delete(setdiff(get(ud.tabs.handles2.hp,'children'),ud.tabs.handles3))
DCM = ud.DCM;
Ep = ud.Ep;
DCM = ud.DCM;
%-graph------------------------------------------------------------
x  = [1:length(DCM.y)]*DCM.Y.dt;
ha = subplot(2,1,1,...
    'parent',ud.tabs.handles2.hp,...
    'nextplot','add');
pos = get(ha,'position');
set(ha,'position',[0.2 pos(2) 0.6 pos(4)])
plot(ha,x,DCM.y(:,i) + DCM.R(:,i),':');
if ud.estimated
    plot(ha,x,DCM.y(:,i));
    title(ha,[DCM.Y.name{i},': response and prediction over time' ],...
        'FontSize',12);
    legend(ha,{'observed','predicted'})
else
    title(ha,['ROI ',DCM.Y.name{i},': response over time' ],...
        'FontSize',12);
end
set(ha,'xlim',[x(1),x(end)])
grid(ha,'on')
box(ha,'off')
ylabel(ha,'signal change (A.U.)')
xlabel(ha,'time (seconds)');

if ud.estimated
    ha = subplot(2,1,2,...
        'parent',ud.tabs.handles2.hp,...
        'nextplot','add');
    pos = get(ha,'position');
    set(ha,'position',[0.2 pos(2) 0.6 pos(4)])
    mi = min([DCM.y(:,i);DCM.y(:,i) + DCM.R(:,i)]);
    ma = max([DCM.y(:,i);DCM.y(:,i) + DCM.R(:,i)]);
    plot(ha,[mi;ma],[mi;ma],'r')
    plot(ha,DCM.y(:,i),DCM.y(:,i) + DCM.R(:,i),'k.');
    set(ha,'xlim',[mi ma],'ylim',[mi ma])
    title(ha,['ROI ',DCM.Y.name{i},': response versus prediction' ],...
        'FontSize',12);
    grid(ha,'on')
    box(ha,'off')
    ylabel(ha,'observed')
    xlabel(ha,'predicted');
    str = ['Variance explained: ',num2str(100*ud.R2.reg(i),'%4.2f'),'%'];
    d = ma-mi;
    ht = text(mi+0.05*d,ma-0.05*d,str,'parent',ha,'color','r','FontSize',12);
end

try;VBA_getSubplots; end



function myA(hfp)
try
    hf = hfp;
catch
    hf = get(get(gco,'parent'),'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles2.hp,'children'))
DCM = ud.DCM;
Ep = ud.Ep;
a(:,:,1) = ud.Pp.A;%DCM.Pp.A;
a(:,:,2) = Ep.A;
ha = displayMat(a,ud.tabs.handles2.hp,DCM,'A',ud.dim);



function myB(hfp)
try
    hf = hfp;
catch
    hf = get(get(gco,'parent'),'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles2.hp,'children'))
DCM = ud.DCM;
if size(DCM.U.u,1) >= 1
    ud.tabs.handles3(1) = uicontrol(...
        'style','popupmenu',...
        'parent',ud.tabs.handles2.hp,...
        'units','normalized',...
        'fontsize',12,...
        'string',DCM.U.name,...
        'callback',@myBi);
    set(ud.tabs.handles3(1),'position',[0.85 0.9 0.12 0.02])
    ud.tabs.handles3(2) = uicontrol(...
        'style','text',...
        'parent',ud.tabs.handles2.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.82 0.93 0.18 0.02],...
        'fontsize',12,...
        'string','display input...');
    set(hf,'userdata',ud)
    feval(@myBi,ud.tabs.handles3(1),[])
else
    ht = uicontrol(...
        'style','text',...
        'parent',ud.tabs.handles2.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.25 0.5 0.5 0.1],...
        'fontsize',12,...
        'string','no input modulatory effect');
end



function myBi(hObject,evt)
hf = get(get(get(hObject,'parent'),'parent'),'parent');
ind = get(hObject,'Value');
ud = get(hf,'userdata');
delete(setdiff(get(ud.tabs.handles2.hp,'children'),ud.tabs.handles3))
DCM = ud.DCM;
Ep = ud.Ep;
if isequal(unique(DCM.b(:,:,ind)),0)
    ht = uicontrol(...
    'style','text',...
        'parent',ud.tabs.handles2.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.25 0.5 0.5 0.1],...
        'fontsize',12,...
        'string',...
        ['no modulatory effect for input ''',DCM.U.name{ind},'''']);
else
    a(:,:,1) = ud.Pp.B(:,:,ind);%DCM.Pp.B(:,:,ind);
    a(:,:,2) = Ep.B(:,:,ind);
    ha = displayMat(a,ud.tabs.handles2.hp,DCM,'B',ud.dim);
end




function myC(hfp)
try
    hf = hfp;
catch
    hf = get(get(gco,'parent'),'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles2.hp,'children'))
DCM = ud.DCM;
Ep = ud.Ep;
if isempty(DCM.c) || isequal(unique(DCM.c(:)),0)
    ht = uicontrol(...
        'style','text',...
        'parent',ud.tabs.handles2.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.25 0.5 0.5 0.1],...
        'fontsize',12,...
        'string','no gain on driving inputs.');
else
    a(:,:,1) = ud.Pp.C;%DCM.Pp.C;
    a(:,:,2) = Ep.C;
    ha = displayMat(a,ud.tabs.handles2.hp,DCM,'C',ud.dim);
end


function myD(hfp)
try
    hf = hfp;
catch
    hf = get(get(gco,'parent'),'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles2.hp,'children'))
DCM = ud.DCM;
ud.tabs.handles3(1) = uicontrol(...
    'style','popupmenu',...
    'parent',ud.tabs.handles2.hp,...
    'units','normalized',...
    'fontsize',12,...
    'string',DCM.Y.name,...
    'callback',@myDi);
set(ud.tabs.handles3(1),'position',[0.85 0.9 0.12 0.02])
ud.tabs.handles3(2) = uicontrol(...
    'style','text',...
    'parent',ud.tabs.handles2.hp,...
    'BackgroundColor',get(hf,'color'),...
    'units','normalized',...
    'position',[0.82 0.93 0.18 0.02],...
    'fontsize',12,...
    'string','display gating...');
set(hf,'userdata',ud)
feval(@myDi,ud.tabs.handles3(1),[])



function myDi(hObject,evt)
hf = get(get(get(hObject,'parent'),'parent'),'parent');
ind = get(hObject,'Value');
ud = get(hf,'userdata');
delete(setdiff(get(ud.tabs.handles2.hp,'children'),ud.tabs.handles3))
DCM = ud.DCM;
Ep = ud.Ep;
if isempty(DCM.d) || isequal(unique(DCM.d(:,:,ind)),0)
    ht = uicontrol(...
        'style','text',...
        'parent',ud.tabs.handles2.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.25 0.5 0.5 0.1],...
        'fontsize',12,...
        'string',...
        ['no nonlinear gating effect for region ''',DCM.Y.name{ind},'''']);
else
    a(:,:,1) = ud.Pp.D(:,:,in);%DCM.Pp.D(:,:,ind);
    a(:,:,2) = Ep.D(:,:,ind);
    ha = displayMat(a,ud.tabs.handles2.hp,DCM,'D',ud.dim);
end


function myCov(hfp)
try
    hf = hfp;
catch
    hf = get(get(gco,'parent'),'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles2.hp,'children'))
% get ticks and tick labels
tick = [0];
ltick = [];
ticklabel = {'A','B','C'};
% A
ltick = [ltick,numel(ud.Ep.A)/2];
tick = [tick,numel(ud.Ep.A)];
% B
ltick = [ltick,tick(end)+numel(ud.Ep.B)/2];
tick = [tick,tick(end)+numel(ud.Ep.B)];
% C
ltick = [ltick,tick(end)+numel(ud.Ep.C)/2];
tick = [tick,tick(end)+numel(ud.Ep.C)];
% D
if ~isempty(ud.Ep.D)
    ltick = [ltick,tick(end)+numel(ud.Ep.D)/2];
    tick = [tick,tick(end)+numel(ud.Ep.D)];
    ticklabel{end+1} = 'D';
end
% hemodynamic params
ltick = [ltick,(tick(end)+(size(ud.Pcorr,1)-tick(end))/2)];
tick = [tick,tick(end)+(size(ud.Pcorr,1)-tick(end))];
ticklabel{end+1} = 'hemo';
tick = tick +0.5;
tick = tick(2:end-1);
ltick = ltick + 0.5;
% display posterior correlation matrix
ha = axes('parent',ud.tabs.handles2.hp);
ud.Pcorr = ud.Pcorr + diag(NaN.*diag(ud.Pcorr));
imagesc(ud.Pcorr,'parent',ha)
try
    set(ha,...
        'xtick',ltick,...
        'ytick',ltick,...
        'xticklabel',ticklabel,...
        'yticklabel',ticklabel,...
        'box','off',...
        'nextplot','add');
    for i=1:length(tick)
        plot(ha,...
            [0.5 size(ud.Pcorr,1)+0.5],...
            [tick(i) tick(i)],...
            'color',[1 1 1])
        plot(ha,...
            [tick(i) tick(i)],...
            [0.5 size(ud.Pcorr,1)+0.5],...
            'color',[1 1 1])
    end
end
grid(ha,'off')
axis(ha,'square')
title(ha,'Parameters posterior correlation matrix','FontSize',12)
set(ha,'clim',[-34/32 1]);
col = colormap('jet');
col(1,:) = 0.5*ones(1,3);
colormap(ha,col);
cb = colorbar('peer',ha);
try;VBA_getSubplots; end



function myResiduals(hfp)
try
    hf = hfp;
catch
    hf = get(get(gco,'parent'),'parent');
end
ud = get(hf,'userdata');
delete(get(ud.tabs.handles2.hp,'children'))
CR = VBA_spm_autocorr(ud.DCM.R');
[nt,nreg] = size(ud.DCM.Y.y);
TR = ud.DCM.Y.dt;
str = cell(nreg,1);
for i=1:nreg
    str{i} = ['ROI #',num2str(i)];
end
ha = subplot(2,1,1,'parent',ud.tabs.handles2.hp);
plot(ha,[TR:TR:nt*TR],ud.DCM.R)
axis(ha,'tight')
title(ha,'residuals time series','fontsize',11)
xlabel(ha,'time (secs)','fontsize',8)
ylabel(ha,'residuals: e(t) = y(y) - g(x(t))','fontsize',8)
legend(ha,str)
set(ha,'ygrid','on','xgrid','off','box','on');

ha = subplot(2,1,2,'parent',ud.tabs.handles2.hp);
plot(ha,[-nt*TR:TR:nt*TR-1],fftshift(CR)')
axis(ha,'tight')
title(ha,'residuals empirical autocorrelation','fontsize',11)
xlabel(ha,'lag tau (secs)','fontsize',8)
ylabel(ha,'Corr[ e(t), e(t+tau) ]','fontsize',8)
legend(ha,str)
set(ha,'ygrid','on','xgrid','off','box','on');
try;VBA_getSubplots; end


function ha = displayMat(a,hParent,DCM,matType,dim)
ha = subplot(2,1,1,'parent',hParent);
imagesc(a(:,:,2),'parent',ha)
title(ha,'posterior mean: E[theta|y,m]','FontSize',12)
if isequal(matType,'C')
    set(ha,...
        'XTick',[1:dim.m],'XTickLabel',DCM.U.name,...
        'YTick',[1:dim.l],'YTickLabel',DCM.Y.name)
    xlabel(ha,'inputs')
    ylabel(ha,'to')
else
    set(ha,...
        'XTick',[1:dim.l],'XTickLabel',DCM.Y.name,...
        'YTick',[1:dim.l],'YTickLabel',DCM.Y.name)
    xlabel(ha,'from')
    ylabel(ha,'to')
end
axis(ha,'square')
dval = max(max(a(:,:,2))) - min(min(a(:,:,2)));
if isequal(dval,0)
    dval = max(max(a(:,:,2)));
    if isequal(dval,0)
        dval = 1e-3;
    end
else
    dval = dval + 1e-4;
end
clim = [min(min(a(:,:,2)))-2*dval/63,max(max(a(:,:,2)))];
set(ha,'clim',sort(clim));
try
    cb = colorbar('peer',ha);
    set(cb,'location','SouthOutside')
catch
    colorbar
end
ha(2) = subplot(2,1,2,'parent',hParent);
imagesc(a(:,:,1),'parent',ha(2))
title(ha(2),['Savage-Dickey ratios: 1-P(theta=0|y,m)'],'FontSize',12)
if isequal(matType,'C')
    set(ha(2),...
        'XTick',[1:dim.m],'XTickLabel',DCM.U.name,...
        'YTick',[1:dim.l],'YTickLabel',DCM.Y.name)
    xlabel(ha(2),'inputs')
    ylabel(ha(2),'to')
else
    set(ha(2),...
        'XTick',[1:dim.l],'XTickLabel',DCM.Y.name,...
        'YTick',[1:dim.l],'YTickLabel',DCM.Y.name)
    xlabel(ha(2),'from')
    ylabel(ha(2),'to')
end
axis(ha(2),'square')
set(ha(2),'clim',[-1/32 1]);
col = colormap('jet');
col(1,:) = 0.5*ones(1,3);
colormap(ha(2),col);
try
    cb = colorbar('peer',ha(2));
    set(cb,'location','SouthOutside')
catch
    colorbar
end
for i = 1:size(a,1)
    for j = 1:size(a,2)
        text(j,i,num2str(a(i,j,2),'%5.2f'),...
            'Parent',ha(1),...
            'FontSize',10,...
            'HorizontalAlignment','Center')
        text(j,i,num2str(a(i,j,1),'%5.2f'),...
            'Parent',ha(2),...
            'FontSize',10,...
            'HorizontalAlignment','Center')
    end
end
try;VBA_getSubplots; end


function R = VBA_spm_autocorr(y)
% computes sample autocorrelation function of signal y
[n,t] = size(y);
R = zeros(n,t*2);
% standardize y
my = mean(y,2);
sy = std(y,[],2);
y = y - repmat(my,1,t);
y = diag(1./sy)*y;
% compute auto-correlation using FFT
for i=1:n
    fr = fft(y(i,:),t*2);
    tfr = fr';
    S = fr(:).*tfr(:);
    R(i,:) = ifft(S)'/t;
end



function DCM = myspm_dcm_test(DCM)
% tests whether DCM params are zero using Savage-Dickey ratios
n = size(DCM.Ep.A,1);
nu = size(DCM.Ep.C,2);
E = VBA_spm_vec(DCM.Ep);
V = DCM.Cp;
E0 = VBA_spm_vec(DCM.M.pE);
V0 = DCM.M.pC;
sdr.A = zeros(n,n);
for i=1:n
    for j=1:n
        tmp = DCM.Vp;
        tmp.A(i,j) = 0;
        ind = find(full(VBA_spm_vec(tmp))~=full(VBA_spm_vec(DCM.Vp)));
        E0r = E0;
        E0r(ind) = 0;
        V0r = V0;
        V0r(ind,:) = 0;
        V0r(:,ind) = 0;
        if isempty(ind) % this param was fixed
            if isequal(E0,E0r)
                dF = Inf;
            else
                dF = -Inf;
            end
        else
            dF = VBA_spm_log_evidence(E,V,E0,V0,E0r,V0r);
            if isempty(dF)
                dF = NaN;
            end
        end 
        sdr.A(i,j) = 1./(1+exp(dF));
    end
end
sdr.B = zeros(n,n,nu);
for i=1:n
    for j=1:n
        for k=1:nu
            tmp = DCM.Vp;
            tmp.B(i,j,k) = 0;
            ind = find(full(VBA_spm_vec(tmp))~=full(VBA_spm_vec(DCM.Vp)));
            E0r = E0;
            E0r(ind) = 0;
            V0r = V0;
            V0r(ind,:) = 0;
            V0r(:,ind) = 0;
            if isempty(ind) % this param was fixed
                if isequal(E0,E0r)
                    dF = Inf;
                else
                    dF = -Inf;
                end
            else
                dF = VBA_spm_log_evidence(E,V,E0,V0,E0r,V0r);
                if isempty(dF)
                    dF = NaN;
                end
            end
            sdr.B(i,j,k) = 1./(1+exp(dF));
        end
    end
end
sdr.C = zeros(n,nu);
for i=1:n
    for j=1:nu
        tmp = DCM.Vp;
        tmp.C(i,j) = 0;
        ind = find(full(VBA_spm_vec(tmp))~=full(VBA_spm_vec(DCM.Vp)));
        E0r = E0;
        E0r(ind) = 0;
        V0r = V0;
        V0r(ind,:) = 0;
        V0r(:,ind) = 0;
        if isempty(ind) % this param was fixed
            if isequal(E0,E0r)
                dF = Inf;
            else
                dF = -Inf;
            end
        else
            dF = VBA_spm_log_evidence(E,V,E0,V0,E0r,V0r);
            if isempty(dF)
                dF = NaN;
            end
        end
        sdr.C(i,j) = 1./(1+exp(dF));
    end
end
nD = size(DCM.Ep.D,3);
if nD > 0
    for i=1:n
        for j=1:n
            for k=1:n
                tmp = DCM.Vp;
                tmp.D(i,j,k) = 0;
                ind = find(full(VBA_spm_vec(tmp))~=full(VBA_spm_vec(DCM.Vp)));
                E0r = E0;
                E0r(ind) = 0;
                V0r = V0;
                V0r(ind,:) = 0;
                V0r(:,ind) = 0;
                if isempty(ind) % this param was fixed
                    if isequal(E0,E0r)
                        dF = Inf;
                    else
                        dF = -Inf;
                    end
                else
                    dF = VBA_spm_log_evidence(E,V,E0,V0,E0r,V0r);
                    if isempty(dF)
                        dF = NaN;
                    end
                end
                sdr.D(i,j,k) = 1./(1+exp(dF));
            end
        end
    end
end
DCM.sdr = sdr;




function [h,g] = myspm_dcm_graph(xY,A,ha)
% Region and anatomical graph display
% see spm_dcm_graph.m
try; ha; catch; ha = []; end
col   = {'b','g','r','c','m','y','k','w'};
m     = size(xY,2);
L     = [];
S     = [];
for i = 1:m
    L       = [L xY(i).xyz];
    name{i} = xY(i).name(1:min(end,3));
    S       = [S xY(i).spec];
end
if isempty(ha)
    ha(1) = subplot(2,1,1);
    set(ha(1),'position',[0 .5 1 .5])
    ha(2) = 0;
end
cla(ha(1));
% meshsurf = fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii');
meshsurf = fullfile(spm('Dir'),'canonical','cortex_8196.surf.gii');
H = VBA_spm_mesh_render('Disp',meshsurf,struct('parent',ha(1)));
options.query = 'dummy';
options.hfig  = H.figure;
options.ParentAxes = H.axis;
options.handles.ParentAxes = H.axis;
options.markersize = 32;
h = VBA_spm_eeg_displayECD(L,[],8,name,options);
grid(H.axis,'on')
h.handles.hfig = H.figure;
h.handles.ParentAxes  = H.axis;
h.handles.mesh = H.patch;
for i = 1:m
    set(h.handles.ht(i),'FontWeight','bold')
end
set(h.handles.mesh,'FaceAlpha',0.20);
hh = get(intersect(get(get(H.patch,'uiContextMenu'),'children'),findobj('Label','Transparency')),'children');
set(hh(5),'Checked','Off')
set(hh(1),'Checked','On')
W     = max(abs(A),abs(A'));
W     = W - diag(diag(W));
W     = 3*W/max(W(:));
W     = W.*(W > 1/128);
for i = 1:length(A)
    for j = (i + 1):length(A)
        if W(i,j)
            if abs(A(i,j)) > abs(A(j,i)), c = j; else, c = i; end
            line(L(1,[i j]),L(2,[i j]),L(3,[i j]),...
                'Color',col{c},...
                'LineStyle','-',...
                'LineWidth',W(i,j),...
                'Parent',ha(1));
        end
    end
end


function DCM = fillinDcm(DCM)
% fill in missing entries in DCM structure for review
DCM.M.m = size(DCM.c,2); % number of inputs
DCM.M.n = 5*size(DCM.Y.y,2); % number of hidden states
DCM.M.l = size(DCM.Y.y,2); % number of regions
DCM.y = zeros(size(DCM.Y.y)); % predicted data
DCM.R = DCM.Y.y; % residuals of the model
DCM.F = 'this model has not been estimated yet.';
DCM.Cp = zeros(DCM.M.l^2*(1+DCM.M.m+DCM.M.l)+DCM.M.l*(DCM.M.m+5));
DCM.Ep.A = DCM.a;
DCM.Ep.B = DCM.b;
DCM.Ep.C = DCM.c;
DCM.Ep.D = DCM.d;
DCM.Pp.A = 0.5*ones(size(DCM.a));
DCM.Pp.B = 0.5*ones(size(DCM.b));
DCM.Pp.C = 0.5*ones(size(DCM.c));
DCM.Pp.D = 0.5*ones(size(DCM.d));
