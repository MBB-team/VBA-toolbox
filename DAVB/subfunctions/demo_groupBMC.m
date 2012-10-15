% this demo checks the behaviour group BMC

clear all
close all


K = 14; % # models
n = 32; % # subjects

N = 128;
ep1 = zeros(K,N);
ep2 = zeros(K,N);


hf = figure;
ha = axes('parent',hf,'nextplot','add');

for ii=1:N
    
    ii
    
    bias = [2*rand*randn(floor(K/2),1);zeros(ceil(K/2),1)];
    for i=1:n
        L(:,i) = 1*randn(K,1) + bias;
    end
    options.TolFun = 1e-4;
    options.DisplayWin = 0;
    [posterior,out] = VBA_groupBMC(L,options);
    % VBA_displayGroupBMC(posterior,out)
    
    alpha0 = out.options.priors.a';
    [exp_r,xp,r_samp,g_post] = spm_BMS_gibbs(L',alpha0,1e4);
    
    ep1(:,ii) = out.ep;
    ep2(:,ii) = xp;
    
    plot(ha,ep1(:,ii),ep2(:,ii),'k.');
    drawnow
    pause(1)
    
end

return

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



% refit the estimated freqs with probit model for exceedance probs
g_fname = @g_odds2;
dim.n_phi = K;
dim.n = 0;
dim.n_theta = 0;
y = out.Ef;
opt.priors.iQy{1} = pinv(out.Vf);
opt.updateHP = 0;
[p2,o2] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,opt);
ep = VBA_ExceedanceProb(p2.muPhi,p2.SigmaPhi);




