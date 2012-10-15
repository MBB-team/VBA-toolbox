% demo for the Kilner heuristic

% close all
clear all

% compute dispersive propagator
gt = 1e-2:1e-3:2e-1;
gr = -4e-2:1e-3:4e-2;
n1 = length(gt);
n2 = length(gr);
G = zeros(n1,n2);
P = [3e-1;3e-3;2e3];
% for i=1:n2
%     [G(:,i)] = dispersivePropagator(gr(i),gt(:),P);
% end
% hf = figure('color',[1 1 1 ],'name','dispersive propagator');
% ha = subplot(2,1,1,'parent',hf);
% imagesc(G,'parent',ha)
% colormap(flipud(autumn))
% ha = subplot(2,1,2,'parent',hf);
% [X,Y] = meshgrid(gt,gr);
% h = mesh(X,Y,G','parent',ha);
% set(h,'facecolor','flat','edgecolor',0.8*[1 1 1],'edgealpha',0.5)
% colormap(flipud(autumn))
% i0 = find(gr==0);
% i1 = find(gt==2e-2);
% hold on
% plot3(gt,0.*gt,G(:,i0)','color',[1 1 1])
% plot3(gt(i1)*ones(size(gr)),gr,G(i1,:)','color',[1 1 1])
% plot3(gt,5e-2*ones(size(gt)),G(:,i0)','k--')
% hp = plot3(-2e-2.*ones(size(gr)),gr,G(1,:),'k--');
% set(ha,'ytick',[-4e-2:2e-2:4e-2],'yticklabel',[-4:2:4])
% set(ha,'xtick',[0:5e-2:2e-1],'xticklabel',[0:50:200])
% set(ha,'ylim',[-4e-2,5e-2],'xlim',[-2e-2,2e-1])
% xlabel('time (msec)')
% ylabel('distance (cm)')
% zlabel('density of connections')

% excitatory maximum post-synaptic depolarization
P.me = 1e-0*8;
% inhibitory maximum post-synaptic depolarization
P.mi = 1e-0*32;
% excitatory post-synaptic time constant
P.Ke = 1e3/4;
% inhibitory post-synaptic time constant
P.Ki = 1e3/28;
% amplitude of intrinsic connectivity kernels
P.a = 1e3.* [ 0 0 2
              0 0 8
              2 1 0 ];
% intrinsic connectivity decay constant
P.c = 1e1*0.32.*ones(3,3); 
% conduction velocity
P.v = 3.*ones(3,3);
% radius of cortical source
P.l = 50.*1e-3;
% parameters of the Gaussian observation filter
P.phi = [1;1*P.l];
% parameters of the Gaussian observation filter
P.sig = struct('r',0.54,'eta',30*1e-0,'g',0.135);
P.i1 = 3;
P.i2 = 1;
% frequency grid
gridw = 2.^[-2:0.1:10];%.1:1e-1:120;

gridx = 30;%0:1e0:50;
nt = numel(gridx);
gy = zeros(numel(gridw),nt);
ds = zeros(nt,1);
for t=1:nt
    [ds(t)] = dsdv(gridx(t),P);
    [gy(:,t)] = spectralPower3(P,gridw,ds(t))';
%     t
end

% frequency bands
ifb{1} = find(gridw>=0.1 & gridw <=4); % delta 
ifb{2} = find(gridw>=4 & gridw <=8); % theta 
ifb{3} = find(gridw>=8 & gridw <=12); % alpha 
ifb{4} = find(gridw>=12 & gridw <=30); % beta 
ifb{5} = find(gridw>=30 & gridw <=100); % gamma


ngy = abs(gy).^2;
mf = zeros(1,nt);
sp = zeros(1,nt);
pfb = zeros(5,nt);
for t=1:nt
    sp(t) = sum(ngy(:,t));
    ngy(:,t) = ngy(:,t)./sum(ngy(:,t));
    for i=1:5
        pfb(i,t) = sum(ngy(ifb{i},t));
    end
    mf(t) = sum(sqrt(ngy(:,t)).*gridw')/(2*pi);
end
for i=1:5
    pfb(i,:) = pfb(i,:)./sum(pfb(i,:));
end

hf = figure('color',[1 1 1],'name','EEG frequency modulation');

ha = subplot(3,2,1,'parent',hf);
plot(ha,gridx,ds)
title(ha,'ds(z)/dz')
xlabel(ha,'fundamental mode dV (mV)')

ha = subplot(3,2,2,'parent',hf);
plot(ha,gridx,mf);
xlabel(ha,'fundamental mode dV (mV)')
title(ha,'centre frequency (Hz)')


ha = subplot(3,2,3,'parent',hf);
hm = mesh(ha,log(abs(gy).^2));
set(hm,'facecolor','flat','edgecolor',0.8*[1 1 1],'edgealpha',0.5)
axis(ha,'tight')
inf = find(gridw==round(gridw));
% set(ha,'xticklabel',gridx(get(ha,'xtick')))
set(ha,'ytick',inf,'yticklabel',gridw(inf))
xlabel(ha,'fundamental mode dV (mV)')
ylabel(ha,'frequency (rad/s)')
title('log frequency power')

ha = subplot(3,2,4,'parent',hf);
hm = mesh(ha,ngy);
set(hm,'facecolor','flat','edgecolor',0.8*[1 1 1],'edgealpha',0.5)
axis(ha,'tight')
inf = find(gridw==round(gridw));
% set(ha,'xticklabel',gridx(get(ha,'xtick')))
set(ha,'ytick',inf,'yticklabel',gridw(inf))
xlabel(ha,'fundamental mode dV (mV)')
ylabel(ha,'frequency (rad/s)')
title('normalized frequency power')
colormap(flipud(autumn))

ha = subplot(3,2,5,'parent',hf);
plot(ha,gridx,sp);
xlabel(ha,'fundamental mode dV (mV)')
title(ha,'mean frequency power (A.U.)')

ha = subplot(3,2,6,'parent',hf);
plot(ha,gridx,pfb');
xlabel(ha,'fundamental mode dV (mV)')
title(ha,'band frequency power (A.U.)')
legend({'delta','theta','alpha','beta','gamma'})

getSubplots


f_fname = @f_modek;
g_fname = @g_Id;
n_t = 1e2;
u = zeros(1,n_t);
u(2) = 1;
x0 = zeros(3,1);
options.inF = P;
options.inF.dt = 1e-3;
options.inF.C = [0;0;1];
options.inF.k = 1;
options.inF.z_10 = 0;
options.dim = struct('n',3,'n_phi',0,'n_theta',0,'n_t',n_t,'p',3,'n_u',1);
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,[],[],u,Inf,Inf,options,x0);

displaySimulations(y,x,eta,e)





