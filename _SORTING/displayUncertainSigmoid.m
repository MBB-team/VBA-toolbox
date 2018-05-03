function [handles] = displayUncertainSigmoid(posterior,out,hf)
% display analysis results

try
    handles.hf0 = hf;
catch
    handles.hf0 = figure('color',ones(1,3));
end
handles.ha0 = subplot(2,1,1,'parent',handles.hf0,'nextplot','add');
title(handles.ha0,'experimental response')
xlabel(handles.ha0,'u: design control variable (motion coherency)')
ylabel(handles.ha0,'detection probability')
grid(handles.ha0,'on')

y = out.y;
u = out.u;


% check sigmoid-like response curve
ne = 8; % number of bins
mi = min(u);
ma = max(u);
dx = (ma-mi)./ne;
edges = mi:dx:ma;
ne = length(edges)-1;
gx = zeros(ne,1);
stdg = zeros(ne,1);
gridx = zeros(ne,1);
for i=1:ne
    ind = find(u>=edges(i)&u<=edges(i+1));
    gx(i) = mean(y(ind));
    stdg(i) = std(y(ind));
    gridx(i) = mean(edges([i,i+1]));
end
errorbar(gridx,gx,stdg,'k.','parent',handles.ha0)

% add model fit
g0 = 1/2;
slope = exp(posterior.muPhi(1))./4;
vslope = slope.^2.*posterior.SigmaPhi(1,1);
options = out.options;
options.priors = posterior;
[gx,vy] = VBA_getLaplace(sort(u),[],options.g_fname,out.dim,options,0,'diag');
plotUncertainTimeSeries(gx(:)',vy(:)',sort(u)',handles.ha0);
plot(handles.ha0,sort(u),gx,'b.')
yy = [0 g0 1];
xx = (yy-g0)./slope;
vyy = ((xx.*(slope+sqrt(vslope))+g0)-yy).^2;
[haf,handles.hf,handles.hp] = plotUncertainTimeSeries(yy,vyy,xx+posterior.muPhi(2),handles.ha0);
set(handles.hf,'facecolor',[1 0 0])
set(handles.hp,'color',[1 0 0])
plot(handles.ha0,posterior.muPhi(2),g0,'go')
sip = sqrt(posterior.SigmaPhi(2,2));
plot(handles.ha0,[posterior.muPhi(2)-sip,posterior.muPhi(2)+sip],[g0 g0],'g')
legend(handles.ha0,...
    {'binned responses',...
    'sigmoid estimate',...
    '1 sigmoid std',...
    'data samples',...
    'sigmoid slope estimate',...
    '1 sigmoid slope std',...
    'inflexion point estimate',...
    '1 inflexion point std'})
box(handles.ha0,'off')
set(handles.ha0,'xgrid','off')

[ny,nx] = hist(u);
handles.ha02 = subplot(2,1,2,'parent',handles.hf0);
bar(handles.ha02,nx,ny,'facecolor',[.8 .8 .8]);
xl0 = get(handles.ha0,'xlim');
set(handles.ha02,'xlim',xl0);
xlabel(handles.ha02,'u: design control variable (stimulus contrast)')
ylabel(handles.ha02,'empirical distribution of u')
title(handles.ha02,'empirical histogram')
grid(handles.ha02,'on')
box(handles.ha02,'off')
set(handles.ha02,'xgrid','off')

try getSubplots; end
