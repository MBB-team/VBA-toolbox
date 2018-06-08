function [] = VBA_spm_dcm_buid(DCM)

% build GUI for DCM results reviewing

if ~nargin
    [P, sts] = VBA_spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    if ~sts, return; end
end
if ~isstruct(P)
    try
        load(P);
        DCM.filename = P;
    catch
        error('Cannot load DCM file "%s".',P);
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
    ud.display.U  = 0.9; % p-value threshold for display
    ud.estimated = 1;
catch
    DCM = fillinDcm(DCM);
    ud.dim.m  = DCM.M.m; % number of inputs
    ud.dim.n  = DCM.M.n; % number of hidden states
    ud.dim.l  = DCM.M.l; % number of regions (responses)
    ud.display.U  = 0.9; % p-value threshold for display
    ud.estimated = 0;
end
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
end
ud.Pcorr = VBA_cov2corr(DCM.Cp);
pos0 = get(0,'screenSize');
pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];
[pathstr,filename,ext] = fileparts(DCM.filename);
figname = ['Explore ',filename];
hfp = figure(...
    'position',pos,...
    'color',[1 1 1],...
    'name',figname,...
    'tag','exploreDCM',...
    'Renderer','OpenGL',...
    'menu','none');
uimenu(hfp,'Label','Open new DCM file','Callback','VBA_spm_dcm_explore');
if ud.estimated
    labels = {'summary','parameters','kernels'};
    callbacks = {@mySummary,@myNeuro,@myHemo};
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
[ud.tabs.handles2] = VBA_spm_uitab(ud.tabs.handles.hp,...
    labels,callbacks,'exploreDCM0',1);
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
    labels = {'A','B','C','D','cov'};
    callbacks = {@myA,@myB,@myC,@myD,@myCov};
else
    labels = {'A','B','C','D'};
    callbacks = {@myA,@myB,@myC,@myD};
end
[ud.tabs.handles2] = VBA_spm_uitab(ud.tabs.handles.hp,...
    labels,callbacks,'exploreDCM0',1);
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
ud = get(hf,'userdata');
ud.tabs.handles3(1) = uicontrol('style','popupmenu',...
    'parent',ud.tabs.handles.hp,...
    'units','normalized',...
    'position',[0.85 0.9 0.12 0.02],...
    'fontsize',12,...
    'string',DCM.U.name,...
    'callback',@myKernel);
ud.tabs.handles3(2) = uicontrol('style','text',...
    'parent',ud.tabs.handles.hp,...
    'BackgroundColor',get(hf,'color'),...
    'units','normalized',...
    'position',[0.82 0.93 0.18 0.02],...
    'fontsize',12,...
    'string','display input...');
set(hf,'userdata',ud)
feval(@myKernel,ud.tabs.handles3(1),[])



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
ha = subplot(2,1,1,...
    'parent',ud.tabs.handles.hp);
plot(ha,x,y)
set(ha,'XLim',[1 max(x)])
pos = get(ha,'position');
set(ha,'position',[0.2 pos(2) 0.6 pos(4)])
title(ha,...
    ['neuronal impulse responses to input ''' DCM.U.name{ind} ''''],...
    'fontsize',12)
grid(ha,'on')
xlabel(ha,'time (seconds)')
legend(ha,DCM.Y.name)

% input effects - hemodynamic
%--------------------------------------------------------------
y = DCM.K1(:,:,ind);
k = DCM.H1(:,:,ind);
ha = subplot(2,1,2,...
    'parent',ud.tabs.handles.hp,...
    'nextplot','add');
pos = get(ha,'position');
set(ha,'position',[0.2 pos(2) 0.6 pos(4)])
plot(ha,x,k)
plot(ha,x,y,':')
set(ha,'XLim',[1 max(x)])
title(ha,...
    ['hemodynamic impulse responses to input ''' DCM.U.name{ind} ''''],...
    'fontsize',12)
grid(ha,'on')
xlabel(ha,'time (seconds)')
legend(ha,DCM.Y.name)
try; VBA_getSubplots (); end


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
str{3} = sprintf(['Dimensions of the model:','\n ',...
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
    str{6} = sprintf(['    - observation function: ',gfn,'\n',...
        '    - evolution function: ',ffn,'\n',...
        '    - integration function: ',ifn,'\n ']);
else
    str{6} = [];
end
na = sum(DCM.a(:)) - sum(diag(DCM.a));
nb = sum(DCM.b(:));
nc = sum(DCM.c(:));
nd = sum(DCM.d(:));
str{7} = sprintf(['DCM structure:','\n ',...
    '    - number of connections: ',num2str(na),...
    ' (+',num2str(DCM.M.l),' inhibitory self-connections)','\n ',...
    '    - number of modulatory effects: ',num2str(nb),'\n ',...
    '    - number of input-state couplings: ',num2str(nc),'\n ',...
    '    - number of (nonlinear) gating effects: ',num2str(nd),'\n ']);

if isnumeric(DCM.F)
    str{8} = ['Log model evidence: log p(y|m) > ',num2str(DCM.F,'%4.3e')];
else
    try
        str{8} = ['Warning: ',DCM.F];
    catch
        str{8} = ' ';
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
    ha = subplot(2,1,1,...
        'parent',ud.tabs.handles2.hp,...
        'visible','off');
    VBA_spm_dcm_display(DCM.xY,[],[],ha);
    colormap(ha,'gray')
    hc = get(ha,'children');
    if isempty(hc)
        delete(ha)
    else
        set(ha,'visible','on','xtick',[],'ytick',[])
    end
    ha = subplot(2,1,2,...
        'parent',ud.tabs.handles2.hp);
    line([0 4],[0 0],'parent',ha)
    text(0,-1,'Name',...
        'FontSize',14,...
        'parent',ha)
    text(2,-1,'Voxels',...
        'FontSize',14,...
        'parent',ha)
    text(3,-1,'Location (mm)',...
        'FontSize',14,...
        'parent',ha)
    line([0 4],[-2 -2],...
        'LineWidth',4,...
        'parent',ha)
    y = -3;
    for i = 1:length(DCM.xY)
        name = DCM.xY(i).name;
        N    = length(DCM.xY(i).s);
        L    = DCM.xY(i).xyz;
        r    = DCM.xY(i).spec;
        text(0,y,name,...
            'FontWeight','bold',...
            'FontSize',12,...
            'parent',ha)
        text(2,y,sprintf('%0.0f',N),...
            'FontSize',12,...
            'parent',ha)
        text(3,y,sprintf('%-4.0f %-4.0f %-4.0f',L),...
            'FontSize',12,...
            'parent',ha)
        y = y - 1;
    end
    line([0 4],[y y],'parent',ha)
    axis(ha,'off')
    try;VBA_getSubplots; end
else
    ht = uicontrol('style','text',...
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
for i = 1:ud.dim.m
    ha = subplot(ud.dim.m,1,i,...
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
    title(ha,['ROI ',DCM.Y.name{i},': response and prediction over time' ],...
        'FontSize',12);
    legend(ha,'predicted', 'observed')
else
    title(ha,['ROI ',DCM.Y.name{i},': response over time' ],...
        'FontSize',12);
end
set(ha,'xlim',[x(1),x(end)])
grid(ha,'on')
box(ha,'on')
ylabel(ha,'signal change (A.U.)')
xlabel(ha,'time (scans)');

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
    box(ha,'on')
    ylabel(ha,'observed')
    xlabel(ha,'predicted');
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
a(:,:,1) = DCM.Pp.A;
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



function myBi(hObject,evt)
hf = get(get(get(hObject,'parent'),'parent'),'parent');
ind = get(hObject,'Value');
ud = get(hf,'userdata');
delete(setdiff(get(ud.tabs.handles2.hp,'children'),ud.tabs.handles3))
DCM = ud.DCM;
Ep = ud.Ep;
if isequal(unique(DCM.b(:,:,ind)),0)
    ht = uicontrol('style','text',...
        'parent',ud.tabs.handles2.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.25 0.5 0.5 0.1],...
        'fontsize',12,...
        'string',...
        ['no modulatory effect for input ''',DCM.U.name{ind},'''']);
else
    a(:,:,1) = DCM.Pp.B(:,:,ind);
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
    ht = uicontrol('style','text',...
        'parent',ud.tabs.handles2.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.25 0.5 0.5 0.1],...
        'fontsize',12,...
        'string','no gain on driving inputs.');
else
    a(:,:,1) = DCM.Pp.C;
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
ud.tabs.handles3(1) = uicontrol('style','popupmenu',...
    'parent',ud.tabs.handles2.hp,...
    'units','normalized',...
    'fontsize',12,...
    'string',DCM.Y.name,...
    'callback',@myDi);
set(ud.tabs.handles3(1),'position',[0.85 0.9 0.12 0.02])
ud.tabs.handles3(2) = uicontrol('style','text',...
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
    ht = uicontrol('style','text',...
        'parent',ud.tabs.handles2.hp,...
        'BackgroundColor',get(hf,'color'),...
        'units','normalized',...
        'position',[0.25 0.5 0.5 0.1],...
        'fontsize',12,...
        'string',...
        ['no nonlinear gating effect for region ''',DCM.Y.name{ind},'''']);
else
    a(:,:,1) = DCM.Pp.D(:,:,ind);
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
title(ha,...
    'Parameters posterior correlation matrix','FontSize',12)
set(ha,'clim',[-34/32 1]);
col = colormap('jet');
col(1,:) = 0.5*ones(1,3);
colormap(ha,col);
cb = colorbar('peer',ha);
try;VBA_getSubplots;end


function ha = displayMat(a,hParent,DCM,matType,dim)
ha = subplot(2,1,1,'parent',hParent);
imagesc(a(:,:,2),'parent',ha)
title(ha,'fixed effects','FontSize',12)
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


ha(2) = subplot(2,1,2,...
    'parent',hParent);
imagesc(a(:,:,1),'parent',ha(2))
title(ha(2),['posterior probabilities: P(|effect|>0)'],'FontSize',12)
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

