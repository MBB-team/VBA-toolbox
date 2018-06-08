% this demo looks at analytical approximations to the product of two
% sigmoids

close all
clear all

x = -50:5e-2:50;
nx = length(x);
s1 = VBA_sigmoid(x);

a = -2:4; % slope
b = 0;%-4:4; % inflexion point
na = length(a);
nb = length(b);

dim.p = nx;
dim.n_t =1;
dim.n_phi = 2;
dim.n_theta = 0;
dim.n = 0;
g_fname = @g_sig_u;

priors.muPhi = zeros(dim.n_phi,1);
priors.SigmaPhi = 1e4*eye(dim.n_phi);
options.priors = priors;
options.DisplayWin = 0;

P = zeros(na,nb,2);
hf = figure('color',[1 1 1]);
ha = axes('parent',hf);

for i=1:na
    for j=1:nb
        s2 = VBA_sigmoid(x, 'slope', exp(a(i)), 'center', b(j));
        s12 = s1.*s2;
        cla(ha)
        plot(ha,x,s1,'b--')
        hold(ha,'on')
        plot(ha,x,s2,'g--')
        plot(ha,x,s12,'r')
        drawnow
        [posterior,out] = VBA_NLStateSpaceModel(VBA_vec(s12),VBA_vec(x),[],g_fname,dim,options);
        P(i,j,1) = posterior.muPhi(1);
        P(i,j,2) = posterior.muPhi(2);
%         pause
    end
end

hf = figure('color',[1 1 1]);
ha = subplot(2,2,1,'parent',hf);
imagesc(P(:,:,1)','parent',ha)
colorbar('peer',ha)
title(ha,'slope')
xlabel(ha,'slope')
ylabel(ha,'inflextion point')
set(ha,'xtick',1:na,'xticklabel',a,'ytick',1:nb,'yticklabel',b)
ha = subplot(2,2,2,'parent',hf);
imagesc(P(:,:,2)','parent',ha)
colorbar('peer',ha)
title(ha,'inflexion point')
xlabel(ha,'slope')
ylabel(ha,'inflextion point')
set(ha,'xtick',1:na,'xticklabel',a,'ytick',1:nb,'yticklabel',b)

