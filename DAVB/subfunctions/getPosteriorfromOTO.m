function [] = getPosteriorfromOTO(posterior,out)
% OTO: recovers posterior sufficient statistics from VB inversion

mux = posterior.muX(1:2,:);
sx = [mux(1,:).*(1-mux(1,:));posterior.muX(3,:)];

V = getVar(posterior.SigmaX.current);
V = [zeros(1,size(V,2));V(2,:)];

plotUncertainTimeSeries(mux,sx);
title('Without accounting for OTO posterior variance')
legend({'outcome identity','log cue-outcome association'})
plotUncertainTimeSeries(mux,sx+V);
title('Accounting for OTO posterior variance')
legend({'outcome identity','log cue-outcome association'})

Var = sx + V;
er1 = mux(2,:) - sqrt(Var(2,:));
er2 = mux(2,:) + sqrt(Var(2,:));
figure
plot(sigm(mux(2,:)));
hold on
yp = [sigm(er2),fliplr(sigm(er1))];
xp = [1:size(mux,2),fliplr(1:size(mux,2))];
hf = fill(xp,yp,'r',...
    'facecolor','b',...
    'edgealpha',0,...
    'facealpha',0.25);
title('cue-outcome association passed through the SIGMOID mapping')

function VX = getVar(SigmaX)
T = length(SigmaX);
n = size(SigmaX{1},1);
VX = zeros(n,T);
for t=1:T
    VX(:,t) = diag(SigmaX{t});
end

