function [out] = addBestLinearPredictor(ho,verbose)
% fits a GLM on current graphical object (and adds the line on the graph)
try,ho;catch,ho=gco;end
try,verbose;catch,verbose=1;end
if ~isequal(get(ho,'type'),'line')
    disp('addbestLinearPredictor: current graphical object is not a line plot!')
    out = [];
    return
end
out.gco = ho;
try, verbose; catch; verbose = 0; end
x = vec(get(ho,'xdata'));
n = length(x);
X = [x,ones(n,1)];
y = vec(get(ho,'ydata'));
[pv,stat,df,out] = GLM_contrast(X,y,[1;0],'F',verbose);
if ~verbose
    out.pv = pv;
    out.stat = stat;
    out.df = df;
end
try
    set(out.handles.hf,'name','addbestLinearPredictor on current graphical object')
end
ha = get(ho,'parent');
status = get(ha,'nextplot');
xlim = get(ha,'xlim');
ylim = get(ha,'ylim');
set(ha,'nextplot','add');
col = get(ho,'color');
yhat = [vec(xlim),ones(2,1)]*out.b;
out.hp = plot(ha,xlim,yhat','color',col)
set(ha,'nextplot',status,'xlim',xlim,'ylim',ylim);



