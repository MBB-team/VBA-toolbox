% demo GLM under multiple covariance components
% This demo exemplifies how to wrap the VBA model inversion scheme to make
% it estimate multipe covariance components within the context of a GLM
% analysis.

% function demo_covComp

Qy{1} = kron(diag([1;0;0]),eye(8));
Qy{2} = kron(diag([0;1;0]),eye(8));
Qy{3} = kron(diag([0;0;1]),eye(8));
my = [-2;-4;0];

Qp{1} = kron(diag([1;0]),eye(2));
Qp{2} = kron(diag([0;1]),eye(2));
mp = [-1;-2];

X = randn(24,4);

theta = 0;
for j=1:length(Qp)
    theta = theta + exp(mp(j))*VBA_sqrtm(Qp{j})*randn(4,1);
end

e = 0;
for i=1:length(Qy)
    e = e + exp(my(i))*VBA_sqrtm(Qy{i})*randn(24,1);
end

y = X*theta + e;

[th,eh,hy,hp] = GLM_covComp(y,X,Qy,Qp,0.*my,0.*mp);

hf = figure;
ha = subplot(3,2,1,'parent',hf,'nextplot','add');
bar(th,'parent',ha)
plot(ha,theta,'go')
title(ha,'GLM coef')
ha = subplot(3,2,2,'parent',hf,'nextplot','add');
plot(ha,theta,th,'.')
grid(ha,'on')
title(ha,'GLM coef')
xlabel(ha,'simulated')
xlabel(ha,'estimated')
ha = subplot(3,2,3,'parent',hf,'nextplot','add');
bar([hy;hp],'parent',ha)
plot(ha,exp([my(:);mp(:)]),'go')
title(ha,'hyperparameters')
ha = subplot(3,2,4,'parent',hf,'nextplot','add');
plot(ha,exp([my(:);mp(:)]),[hy;hp],'.')
grid(ha,'on')
title(ha,'hyperparameters')
xlabel(ha,'simulated')
xlabel(ha,'estimated')
ha = subplot(3,2,5,'parent',hf,'nextplot','add');
bar(eh,'parent',ha)
plot(ha,e,'go')
title(ha,'residuals')
ha = subplot(3,2,6,'parent',hf,'nextplot','add');
plot(ha,e,eh,'.')
grid(ha,'on')
title(ha,'residuals')
xlabel(ha,'simulated')
xlabel(ha,'estimated')



