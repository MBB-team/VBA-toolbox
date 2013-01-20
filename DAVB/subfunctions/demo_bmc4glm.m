% This demo uses Monte-Carlo simulations to validate the Bayesian model
% comparison ability in the context of non-nested models.
% NB: the model evidence is evaluated at the frequentist limit.

N = 1e4; % #simulations
d = 16; % data dimension
r = 4; % parameter dimension

s = 1e-1; % noise variance
s2 = 1e0; % signal power

for i = 1:N
    
    X1 = randn(d,r);
    X2 = randn(d,r);%X1(:,2:end);
    
    b = s2*rand(r,1);
    
    y1 = X1*b + sqrt(s)*randn(d,1);
    y2 = X2*b + sqrt(s)*randn(d,1);
    
    [lev1(i,1)] = lev_GLM(y1,X1);
    [lev1(i,2)] = lev_GLM(y1,X2);
    
    [lev2(i,1)] = lev_GLM(y2,X1);
    [lev2(i,2)] = lev_GLM(y2,X2);
    
    
end


[n1,x1] = hist(lev1(:,1)-lev1(:,2));
[n2,x2] = hist(lev2(:,1)-lev2(:,2));

hf = figure('color',[1 1 1]);
ha = axes('parent',hf,'nextplot','add');
bar(ha,x1,n1,'facecolor','r');
bar(ha,x2,n2,'facecolor','b');
legend(ha,{'true = model 1','true = model 2'})
xlabel(ha,'log p(y|m1) - log(y|m2)')
ylabel(ha,'# simulations')


