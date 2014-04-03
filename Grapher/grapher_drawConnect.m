function hs=grapher_drawConnect(node_A,node_B,weight)
%GRAPHER_DRAWCONNECT   Draw a connection between two DCM nodes
%   Inputs :
%       * two nodes structures (having coordinates specified in their field  |pos|) 
%       * a connection weight
%   Outputs
%       * structure provifing handles to the patches (line and enpoint)
%

%% Constants
WIDTH_MIN = 2;
WIDTH_MAX = 3;
WIDTH_WGT = @(w) 4*sqrt(abs(w));

OFFSET = 6 ;
ENDPOINT_RADIUS = 7 ;

%% 

    vector = node_B.pos-node_A.pos;
    
    % move a bit to avoid feedback overlap
    thet = atan2(vector(2),vector(1));
    thet = thet-(pi/2);
    offset = OFFSET*[cos(thet), sin(thet)];
    
    
    xs = [node_A.pos(1) node_B.pos(1)] + offset(1);
    ys = [node_A.pos(2) node_B.pos(2)] + offset(2);
    zs = -2*eps*[1 1];
    
    cs = [weight weight];
    w = min(max(WIDTH_WGT(weight),WIDTH_MIN),WIDTH_MAX);
    h=patch(xs,ys,zs,weight,'FaceColor','flat','FaceVertexCData',cs','CDataMapping','scaled','EdgeColor','flat','LineWidth',w);%,'LineSmoothing','on'

    
    vector2 = [xs(2) ys(2)] - [xs(1) ys(1)];
    vector2 = vector2/norm(vector2);
    
    th = 0:pi/30:2*pi;
    xunit = ENDPOINT_RADIUS * cos(th) + xs(2) - 45*vector2(1);
    yunit = ENDPOINT_RADIUS * sin(th) + ys(2) - 45*vector2(2);
    zunit = -eps*ones(1,numel(th));
    cs = repmat(weight,1,numel(th));
    h2=patch(xunit,yunit,zunit,cs,'FaceColor','flat','FaceVertexCData',cs','CDataMapping','scaled','EdgeColor','flat','LineWidth',.01);%,,'FaceVertexCData',weight,'LineSmoothing','on'

    hs.line=h;
    hs.end=h2;
end