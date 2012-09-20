% this demo checks the behaviour group BMC

clear all
close all


K = 4; % # models
n = 32; % # subjects

bias = randn(K,1);
bias = [ones(floor(K/2),1);zeros(ceil(K/2),1)];

for i=1:n
    L(:,i) = 1*randn(K,1) + bias;
end
options.TolFun = 1e-4;
[posterior,out] = VBA_groupBMC(L,options);

alpha0 = out.options.priors.a';
[exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (L', alpha0,1e4);

plot(out.options.handles.ha(3),exp_r,'go')
plot(out.options.handles.ha(4),xp,'go')

hf = figure;
subplot(2,2,1)
plot(exp_r(:),out.Ef(:),'k.');
hold on
plot([0 1],[0 1],'r')
title('model frequencies')
subplot(2,2,2)
plot(xp(:),out.ep(:),'k.');
hold on
plot([0 1],[0 1],'r')
title('exceedance probabilities')
subplot(2,2,3)
g = g_post';
plot(g(:),posterior.r(:),'k.');
hold on
plot([0 1],[0 1],'r')
title('model attributions')

% VBA_displayGroupBMC(posterior,out)
