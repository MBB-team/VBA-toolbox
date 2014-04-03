function h=grapher_drawNode(node,color)
%GRAPHER_DRAWNODE   Draw a DCM node
%   Inputs :
%       * a node structure having fileds :
%            - |pos|: coordinates of the node
%            - |lbl|: name of the node
%       * an activity level
%   Outputs
%       * structure provifing handles to the patches
%

%% Constants
RADIUS = 40;
NAME_SIZE = 9;

%%
    th = 0:pi/15:2*pi;
    r=RADIUS;
    xunit = r * cos(th) + node.pos(1);
    yunit = r * sin(th) + node.pos(2);
    zunit = eps*ones(1,numel(th));
    cs = repmat(color,1,numel(th));
    h=patch(xunit,yunit,zunit,cs,'FaceColor','flat','FaceVertexCData',cs','CDataMapping','scaled','EdgeColor',[.6 .6 .6],'LineWidth',3);%,,'FaceVertexCData',weight,'LineSmoothing','on'
    text(node.pos(1),node.pos(2),eps,node.lbl,...
    'HorizontalAlignment','center','Color',[.2 .2 .2],...
    'FontWeight','bold','FontSize',NAME_SIZE,...
    'FontName','Verdana');
end