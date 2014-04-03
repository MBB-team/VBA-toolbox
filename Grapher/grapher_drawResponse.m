function h=grapher_drawResponse(response,color)
    r=40;
    xunit = response.pos(1) + r*[-1 -1 1 1]/2;
    yunit = response.pos(2) + r*[-1 1 1 -1]/2;
    zunit = eps*ones(1,4);
    cs = repmat(color,1,4);
    h=patch(xunit,yunit,zunit,cs,'FaceColor','flat','FaceVertexCData',cs','CDataMapping','scaled','EdgeColor',[.6 .6 .6],'LineWidth','none','LineSmoothing','on');%,,'FaceVertexCData',weight,
    text(response.pos(1),response.pos(2)-r,eps,response.lbl,...
    'HorizontalAlignment','center','Color',[.2 .2 .2],...
    'FontWeight','bold','FontSize',8,...
    'FontName','Verdana');
end