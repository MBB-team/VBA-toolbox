function getSubplots
% adds contextmenu to (subplot) axes for whole window display

ha = findobj('Type','axes');
for i=1:length(ha)
    if ~isequal(get(ha(i),'Tag'),'legend')
        hfig = get(ha(i),'parent');
        % this deals with uipanel, etc:
        while ~isequal(get(hfig,'type'),'figure')
            hfig = get(hfig,'parent');
        end
        hcmenu = uicontextmenu('parent',hfig);
        uimenu(hcmenu,...
            'Label','Export axes to figure',...
            'Callback',@getSubplot);
        set(ha(i),'uicontextmenu',hcmenu);
    end
end


function [ha] = getSubplot(ho,tmp)

ha = gca;
ud.op = get(ha,'parent');
ud.pos = get(ha,'position');
ud.fontsize = get(ha,'fontsize');
ud.ha = ha;
stop = 0;
hp = ha;
while ~stop
    hp = get(hp,'parent');
    if isequal(get(hp,'type'),'figure')
        ud.col = get(hp,'colormap');
        stop = 1;
    end
end

% first deal with potential associated colorbars
ud.hclb = [];
hclb = intersect(findobj('tag','Colorbar'),get(get(ha,'parent'),'children'));
for i=1:length(hclb)
    [haa] = findAxes4Colorbar(hclb(i),hclb);
    if isequal(haa,ha)
        ud.hclb = hclb(i);
        ud.clbFontsize = get(hclb(i),'fontsize');
        break
    end
end

% then get potential legends
ud.leg = [];
hleg = intersect(findobj('tag','legend'),get(get(ha,'parent'),'children'));
for i=1:length(hleg)
    try % matlab 7.1?
        tmp = get(hleg(i),'userdata');
        pos = get(hleg(i),'position');
        haa = tmp.PlotHandle;
        if isequal(haa,ha)
            ud.leg.h = hleg(i);
            ud.leg.pos = pos;
            break
        end
    end
end
    

hf = figure(...
    'visible','off',...
    'color',[1 1 1],...
    'colormap',ud.col);
na = axes('parent',hf);
pos = get(na,'position');
delete(na);

set(ha,...
    'parent',hf,...
    'position',pos,...
    'fontsize',16);

if ~isempty(ud.hclb)
    set(ud.hclb,'parent',hf);
    set(ud.hclb,'fontsize',16);
end
if ~isempty(ud.leg)
    set(ud.leg.h,'parent',hf);
end

hc = get(ha,'children');
nc = length(hc);
for i=1:nc
    if isequal(get(hc(i),'type'),'line')
        ud.linewidth(i) = get(hc(i),'linewidth');
        set(hc(i),'linewidth',2);
        ud.MarkerSize(i) = get(hc(i),'MarkerSize');
        set(hc(i),'MarkerSize',14);
    end
end

ht = get(ha,'title');
ud.titleFontsize = get(ht,'fontsize');
set(ht,'fontsize',16);

hlx = get(ha,'xlabel');
ud.xlbFontsize = get(hlx,'fontsize');
set(hlx,'fontsize',14);

hly = get(ha,'ylabel');
ud.ylbFontsize = get(hly,'fontsize');
set(hly,'fontsize',14);

hlz = get(ha,'zlabel');
ud.zlbFontsize = get(hlz,'fontsize');
set(hlz,'fontsize',14);

set(hf,...
    'visible','on',...
    'userdata',ud,...
    'DeleteFcn',@closeFig)


function closeFig(e1,e2)
ud = get(gcf,'userdata');
try
    hfig = ud.op;
    % this deals with uipanel, etc:
    while ~isequal(get(hfig,'type'),'figure')
        hfig = get(hfig,'parent');
    end
    hcmenu = uicontextmenu('parent',hfig);
    uimenu(hcmenu,...
        'Label','Export axes to figure',...
        'Callback',@getSubplot);
    set(ud.ha,...
        'parent',ud.op,...
        'position',ud.pos,...
        'fontsize',ud.fontsize,...
        'uicontextmenu',hcmenu)
    hc = get(ud.ha,'children');
    nc = length(hc);
    for i=1:nc
        if isequal(get(hc(i),'type'),'line')
            set(hc(i),'linewidth',ud.linewidth(i));
            set(hc(i),'MarkerSize',ud.MarkerSize(i));
        end
    end
    ht = get(ud.ha,'title');
    set(ht,'fontsize',ud.titleFontsize);
    hlx = get(ud.ha,'xlabel');
    set(hlx,'fontsize',ud.xlbFontsize);
    hly = get(ud.ha,'ylabel');
    set(hly,'fontsize',ud.ylbFontsize);
    hlz = get(ud.ha,'zlabel');
    set(hlz,'fontsize',ud.zlbFontsize);
    if ~isempty(ud.hclb)
        set(ud.hclb,'parent',ud.op);
        set(ud.hclb,'fontsize',ud.clbFontsize);
    end
    if ~isempty(ud.leg)
        set(ud.leg.h,...
            'parent',ud.op,...
            'position',ud.leg.pos);
    end
end


function [haa] = findAxes4Colorbar(hclbi,hclb)
pos0 = get(hclbi,'position');
ha = intersect(findobj('type','axes'),...
    get(get(hclbi,'parent'),'children'));
ha = setdiff(ha,hclb);
pos = zeros(length(ha),4);
for i=1:length(ha)
    pos(i,:) = get(ha(i),'position');
end
dx1 = pos0(1) - pos(:,1);
dx2 = pos0(1) - (pos(:,1)+pos(:,3));
dy1 = pos0(2) - pos(:,2);
dy2 = pos0(2) - (pos(:,2)+pos(:,4));
switch get(hclbi,'Location')
    case {'North','NorthOutside'}
        ind = find(dy1>=0);
        [tmp,subind] = min(dx1(ind).^2+dy2(ind).^2);
        haa = ha(ind(subind));
    case {'South','SouthOutside'}
        ind = find(dy1<=0);
        [tmp,subind] = min(dx1(ind).^2+dy2(ind).^2);
        haa = ha(ind(subind));
    case {'East','EastOutside'}
        ind = find(dx1>=0);
        [tmp,subind] = min(dy1(ind).^2+dx2(ind).^2);
        haa = ha(ind(subind));
    case {'West','WestOutside'}
        ind = find(dx1<=0);
        [tmp,subind] = min(dy1(ind).^2+dx2(ind).^2);
        haa = ha(ind(subind));
end
        
        
    



