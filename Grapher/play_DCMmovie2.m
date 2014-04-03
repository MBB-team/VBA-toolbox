function play_DCMmovie2(mov,nodes)

% set(0,'DefaultFigureWindowStyle','undocked')
close all
% figure
%% get structure
 
n_nodes = numel(nodes) ;

%% init mov
fig=figure('OuterPosition',[560 300 500 500],...
            'name','DCM Player',...
            'MenuBar','none','NumberTitle','off',...
            'DockControls','off');
hold on

% connections = [];
for i=1:n_nodes
    for j=1:n_nodes
        if mov(1).pattern(j,i)
        connections{i,j} = grapher_drawConnect(nodes(i),nodes(j),1);
        update_connect( connections{i,j},0);
        else
        connections{i,j} = NaN;
        end
    end
end

state = [];
for i=1:n_nodes
    state(i) = grapher_drawNode(nodes(i),0);  
end

colormap(flipud(cbrewer('div','RdBu',100)));
% colormap(lbmap(100,'RedBlue'))
% load cmap;
% colormap(flipud(cmap));
% colormap(fireice)
% colormap(flipud(othercolor('RdBu11')));
xlim([-250 250])
ylim([-250 250])
caxis([-1 1]);
set(gca,'box','off','xcolor','w','ycolor','w')
set(gcf, 'color', [1 1 1])

drawnow

% ui
n_t = length(mov);

%GUI

sli=uicontrol('Style','slider',...
'Min' ,1,'Max',n_t, ...
'Position',[30,10,445,15], ...
'Value', 1,...
'SliderStep',[1/n_t 10/n_t], ...
'BackgroundColor',[1 1 1],...
'CallBack', @slider_Callback,...
'BusyAction','queue'); 

but=uicontrol('Style','togglebutton',...
'Position',[10,10,16,16],...
'Value',0,...
'Min',0,'Max',1,...
'CallBack', @play,...
'BackgroundColor',[.5 0 0]); 

infos=uicontrol('Style','text',...
'Position',[10,30,300,16],...
'String','',...
'HorizontalAlignment','left',...
'BackgroundColor',[1 1 1]); 

timer=uicontrol('Style','text',...
'Position',[317,30,155,16],...
'String','0 s (t=00000)',...
'HorizontalAlignment','right',...
'BackgroundColor',[1 1 1]); 

setappdata(fig,'mov',mov);
setappdata(fig,'nodes',nodes);
setappdata(fig,'state',state);
setappdata(fig,'connections',connections);
setappdata(fig,'sli',sli);
setappdata(fig,'but',but);
setappdata(fig,'infos',infos);
setappdata(fig,'timer',timer);

% return

%% loop
% get timeseries


% for t=1:numel(mov)
%   update_dcm(state,connections,mov(t));
%   drawnow
% end
  

end

function slider_Callback(hObject, eventdata, handles)
    mov = getappdata(gcf,'mov');
    state = getappdata(gcf,'state');
    connections = getappdata(gcf,'connections');
    t=round(get(hObject,'Value'));
    
    connections=update_dcm(state,connections,mov(t)); 
    setappdata(gcf,'connections',connections);
    update_gui(t)
    
end

function play(hObject, eventdata, handles)
    mov = getappdata(gcf,'mov');
    state = getappdata(gcf,'state');
    connections = getappdata(gcf,'connections');
    sli = getappdata(gcf,'sli');

    but = getappdata(gcf,'but');
    start = round(get(sli,'Value'));
    set(but,'BackgroundColor',[0 .5 0]);
    
for t=start:numel(mov)

    update_dcm(state,connections,mov(t));
    update_gui(t)
    
    if get(but,'Value')==0
        set(but,'BackgroundColor',[.5 0 0]);
        return
    end
    
    drawnow

end
end

% function h=grapher_drawNode(node,color)
%     th = 0:pi/30:2*pi;
%     r=40;
%     xunit = r * cos(th) + node.pos(1);
%     yunit = r * sin(th) + node.pos(2);
%     zunit = eps*ones(numel(th),1);
%     cs = repmat(color,1,numel(th));
%     h=patch(xunit,yunit,zunit,cs,'FaceColor','flat','FaceVertexCData',cs','CDataMapping','scaled','EdgeColor',[.6 .6 .6],'LineWidth',3,'LineSmoothing','on');%,,'FaceVertexCData',weight,
%     text(node.pos(1),node.pos(2),eps,node.lbl,...
%     'HorizontalAlignment','center','Color',[.2 .2 .2],...
%     'FontWeight','bold','FontSize',10,...
%     'FontName','Verdana');
% end    
% 
% function h=grapher_drawConnect(node_A,node_B,weight)
%     vector = node_B.pos-node_A.pos;
%         % move a bit to avoid feedback overlap
% 
%     thet = atan2(vector(2),vector(1));
%     thet = thet-(pi/2);
%     offset = 6*[cos(thet), sin(thet)];
%     
%     xs = [node_A.pos(1) node_B.pos(1)] + offset(1);
%     ys = [node_A.pos(2) node_B.pos(2)] + offset(2);
%     zs = [0 0];
%     
%     cs = [weight weight];
%     h=patch(xs,ys,zs,weight,'FaceColor','flat','FaceVertexCData',cs','CDataMapping','scaled','EdgeColor','flat','LineWidth',4,'LineSmoothing','on');
% end

function update_connect(h,weight)
    set(h.line,'FaceVertexCData',weight*ones(numel(get(h.line,'FaceVertexCData')),1));
    set(h.end,'FaceVertexCData',weight*ones(numel(get(h.end,'FaceVertexCData')),1));
end

function update_node(h,color)
    set(h,'CData',color);
end    

function connections=update_dcm(state,connections,frame)
nodes = getappdata(gcf,'nodes');

n_nodes = length(state);
    
    for i=1:n_nodes
       update_node(state(i),frame.activity(i));
       for j=1:n_nodes
           if frame.pattern(j,i)  
            update_connect(connections{i,j},frame.connectivity(j,i));
           end
       end
    end
end

function update_gui(t)
    sli = getappdata(gcf,'sli');
    infos = getappdata(gcf,'infos');
    timer = getappdata(gcf,'timer');
    mov = getappdata(gcf,'mov');

    set(sli,'Value',t);
    set(timer,'String',sprintf('%.0f s (t=%05d)',mov(t).time,t));
    set(infos,'String',sprintf('Inputs = |%s',sprintf(' %2.2f |',mov(t).input)));
end
    