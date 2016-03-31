function MoveAxisToOrigin(ha)
% move axes TO origin
% function MoveAxisToOrigin(ha)
% IN:
%   - ha: axes' handle

% GET TICKS
X=get(ha,'Xtick');
Y=get(ha,'Ytick');
fs = get(ha,'fontsize');

% GET LABELS
XL=get(ha,'XtickLabel');
YL=get(ha,'YtickLabel');

% GET OFFSETS
Xoff=diff(get(ha,'XLim'))./40;
Yoff=diff(get(ha,'YLim'))./40;

% DRAW AXIS LINEs
hold(ha,'on')
plot(get(ha,'XLim'),[0 0],'k');
plot([0 0],get(ha,'YLim'),'k');

% Plot new ticks  
for i=1:length(X)
    plot(ha,[X(i) X(i)],[0 -Yoff],'-k');
end;
for i=1:length(Y)
   plot(ha,[-Xoff, 0],[Y(i) Y(i)],'-k');
end;

% ADD LABELS
text(X,zeros(size(X))-2.*Yoff,XL,'fontsize',fs,'HorizontalAlignment','center','VerticalAlignment','cap','parent',ha);
text(zeros(size(Y))-2.*Xoff,Y,YL,'fontsize',fs,'HorizontalAlignment','right','VerticalAlignment','middle','parent',ha);

box(ha,'off');
% axis(ha,'square');
axis(ha,'off');

