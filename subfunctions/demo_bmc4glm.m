% This demo demonstrates random-effect group BMS in the context of nested
% models. NB: the model evidence is evaluated at the frequentist limit.

close all
clear variables

N = 32; % # subjects
d = 64; % data dimension
r = 4; % parameter dimension

s = 1e0; % noise variance
s2 = 1e0; % signal power

for i = 1:N
    
    X1 = randn(d,r);
    X2 = X1(:,1:2);%randn(d,r);%X1(:,2:end);
    
    b = sqrt(s2)*rand(r,1);
    
    y1 = X1*b + sqrt(s)*randn(d,1);
    y2 = X2*b(1:2) + sqrt(s)*randn(d,1);
    
    [lev1(i,1)] = lev_GLM(y1,X1);
    [lev1(i,2)] = lev_GLM(y1,X2);
    
    [lev2(i,1)] = lev_GLM(y2,X1);
    [lev2(i,2)] = lev_GLM(y2,X2);
    
    
end

% examplifies classical GLM analysis of y|m2
c = eye(4);
c = c(1:2,:);
[pv,stat,df,all] = GLM_contrast(X1,y2,c','F',1);
set(gcf,'name','classical analysis of y|m2')

% display empirical histogram of log- Bayes factors
[n1,x1] = empiricalHist(lev1(:,1)-lev1(:,2));
[n2,x2] = empiricalHist(lev2(:,1)-lev2(:,2));
hf = figure('color',[1 1 1],'name','Monte-Carlo simulations: distribution of log Bayes factors');
ha = axes('parent',hf,'nextplot','add');
plot(ha,x1,n1,'color','r');
plot(ha,x2,n2,'color','b');
legend(ha,{'true = model 1','true = model 2'})
xlabel(ha,'log p(y|m1) - log(y|m2)')
ylabel(ha,'# simulations')

% perform group-BMS on data generated under model 1
L1 = lev1';
[p1,o1] = VBA_groupBMC(L1);
set(o1.options.handles.hf,'name','group BMS: y|m1')

% perform group-BMS on data generated under model 2
L2 = lev2';
[p2,o2] = VBA_groupBMC(L2);
set(o2.options.handles.hf,'name','group BMS: y|m2')


