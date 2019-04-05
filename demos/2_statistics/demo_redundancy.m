% Redundancy reduction principle
% This demo performs BMC on two nested models, one of which is fully
% redundant. The data is then varied on a pre-defined grid, to assess the
% impact of the complexity penalty term as a function of the data.


clear variables
close all

dim1.n_phi = 1;
dim1.n_theta = 0;
dim1.n=0;
dim2.n_phi = 2;
dim2.n_theta = 0;
dim2.n=0;
opt1.inG.X = 1;
opt1.DisplayWin = 0;
opt1.verbose = 0;
opt2.inG.X = [1,1];
opt2.DisplayWin = 0;
opt2.verbose = 0;

gy = 2.^[-8:8];
ny = length(gy);
for i=1:ny
    i
    y = gy(i);
    [p1,o1] = VBA_NLStateSpaceModel(y,[],[],@g_GLM,dim1,opt1);
    set(gcf,'tag','1')
    [p2,o2] = VBA_NLStateSpaceModel(y,[],[],@g_GLM,dim2,opt2);
    dF(i) = o1.F(end) - o2.F(end);
    m1(i) = p1.muPhi;
    s1(i) = p1.SigmaPhi;
    m2(i) = ones(2,1)'*p2.muPhi;
    s2(i) = ones(2,1)'*p2.SigmaPhi*ones(2,1);
    c = VBA_cov2corr(p2.SigmaPhi);
    c2(i) = c(1,2);
end

hf = figure('color',[1 1 1]);
ha = subplot(2,2,1,'parent',hf);
plot(ha,gy,dF,'marker','.')
set(ha,'xscale','log')
ylabel(ha,'log p(y|m1) - log p(y|m2)')
xlabel(ha,'y')
title(ha,'Bayesian model comparison')
ha = subplot(2,2,2,'parent',hf);
plot(ha,gy,c2,'marker','.')
set(ha,'xscale','log')
ylabel(ha,'Corr[x1,x2|y,m2]')
xlabel(ha,'y')
title(ha,'parameter identifiability under m2')
ha = subplot(2,2,3,'parent',hf);
plotUncertainTimeSeries(m1,s1,gy,ha)
set(ha,'xscale','log')
ylabel(ha,'E[x1|y,m1]')
xlabel(ha,'y')
title(ha,'estimated effect under m1')
ha = subplot(2,2,4,'parent',hf);
plotUncertainTimeSeries(m2,s2,gy,ha)
set(ha,'xscale','log')
ylabel(ha,'E[x1+x2|y,m2]')
xlabel(ha,'y')
title(ha,'estimated effect under m2')


