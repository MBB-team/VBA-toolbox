function [] = getPosteriorfromOTO(posterior,out)
% OTO: recovers posterior sufficient statistics from VB inversion

mux = posterior.muX(1:2,:);
sx = [mux(1,:).*(1-mux(1,:));posterior.muX(3,:)];

V = VBA_getVar(posterior.SigmaX.current);
V = [zeros(1,size(V,2));V(2,:)];

ha = plotUncertainTimeSeries(mux,sx);
title(ha,'Without accounting for OTO posterior variance')
legend(ha,{'outcome identity','log cue-outcome association'})
ha = plotUncertainTimeSeries(mux,sx+V);
title(ha,'Accounting for OTO posterior variance')
legend(ha,{'outcome identity','log cue-outcome association'})

Var = sx + V;
er1 = mux(2,:) - sqrt(Var(2,:));
er2 = mux(2,:) + sqrt(Var(2,:));
hf = figure('color',[1 1 1]);
ha = axes('parent',hf,'nextplot','add');
plot(ha,VBA_sigmoid(mux(2,:)));
yp = [VBA_sigmoid(er2),fliplr(VBA_sigmoid(er1))];
xp = [1:size(mux,2),fliplr(1:size(mux,2))];
fill(xp,yp,'r','facecolor','b','edgealpha',0,'facealpha',0.25,'parent',ha);
title(ha,'cue-outcome association passed through the SIGMOID mapping')


