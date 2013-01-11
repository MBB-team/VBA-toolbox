% demo for the Kilner heuristic

close all
clear all

% compute dispersive propagator
dt = 1e-3;
dr = 2e-4;
gt = 1e-3:dt:8e-2;
gr = 1e-3:dr:1.6e-2;
n1 = length(gt);
n2 = length(gr);
G0 = zeros(n1,n2);
Er = zeros(n1,1);
Ger = zeros(n1,1);
pG = G0;
P = [3e-1;3e-3;2e3];
for i=1:n1
    [G0(i,:)] = 2*pi*gr.*dispersivePropagator(gr(:),gt(i),P)';
    pG(i,:) = G0(i,:)./sum(G0(i,:));
    Er(i) = sum(gr.*pG(i,:));
    Ger(i) = dispersivePropagator(Er(i),gt(i),P);
end
[X,Y] = meshgrid(gt,gr);
dummy=10*log10((G0+eps)/max(G0(:)));
dbG=nan(size(dummy));
dbG(dummy>-30)=G0(dummy>-30);
G = dbG/max(G0(:));



hf = figure('color',[1 1 1 ],'name','dispersive propagator');
ha = subplot(2,1,1,'parent',hf);
% imagesc(flipud(log(G)'),'parent',ha)
imagesc(log(G)','parent',ha)
it = find(gr*5e2==floor(gr*5e2));

colorbar
hold on
v = 2.^(-[1:9]);
% [C,hc] = contour(flipud(log(G)'),log(v(:)));
[C,hc] = contour(log(G)',log(v(:)));
set(hc,'color','k')

Ier = ((Er-min(Er))./max(Er))*length(gr)+1;
hp=plot(1:length(gt),Ier,'color',[1 1 1],'parent',ha)

xlabel('time (msec)')
ylabel('distance (mm)')
title('log-density of connections')
% set(ha,'xticklabel',1e3*gt(get(ha,'xtick')),'ytick',it,'yticklabel',1e3*fliplr(gr(it)))
set(ha,'xticklabel',1e3*gt(get(ha,'xtick')),'ytick',it,'yticklabel',1e3*gr(it),'ydir','reverse')
set(ha,'ydir','normal')

ha = subplot(2,1,2,'parent',hf);
h = mesh(X,Y,G','parent',ha);
set(h,'facecolor','flat','edgecolor',0.8*[1 1 1],'edgealpha',0.5)
colormap(flipud(autumn))
i0 = find(gr==min(abs(gr)));
i1 = find(gt==2e-2);
hold on
[C,hc] = contour3(X,Y,G'+1e-4,v(:),'parent',ha);
set(hc,'EdgeColor','k')
plot3(gt,Er,1e-6+Ger/max(G0(:))','color',[1 1 1])

set(ha,'ytick',[0:2e-3:max(gr)],'yticklabel',[0:2e-3:max(gr)]*1e3)
set(ha,'xtick',[0:1e-2:max(gt)],'xticklabel',[0:1e-2:max(gt)]*1e3)
set(ha,'ylim',[0,max(gr)],'xlim',[0,max(gt)],'zlim',[min(G(:)),max(G(:))])
set(ha,'zscale','log','CameraPosition',[0.251629 0.105342 6.06163])
xlabel('time (msec)')
ylabel('distance (mm)')
zlabel('density of connections')

getSubplots


% return

unit = 1; % if 1: mV, if 1e-3: V

% excitatory maximum post-synaptic depolarization
P.me = unit*8;
% inhibitory maximum post-synaptic depolarization
P.mi = unit*32;
% excitatory post-synaptic time constant
P.Ke = 1e3/4;
% inhibitory post-synaptic time constant
P.Ki = 1e3/32;
% amplitude of intrinsic connectivity kernels (#synapses)
P.a = 1e3.* [ 0 0 2
              0 0 8
              2 1 0 ];
% intrinsic connectivity decay constant (spatial scale of lateral connectivity)
P.c = 1e-2*ones(3,3);%3e-3.*ones(3,3); 
% conduction velocity (along connections)
P.v = 1.*ones(3,3);%3.*ones(3,3);
% radius of cortical source
P.l = 1e-2;%5.*1e-3;
% parameters of the Gaussian observation filter
P.phi = [1;1*P.l];
% parameters of the Gaussian observation filter
P.sig = struct('r',0.54,'eta',30*unit,'g',0.135);
P.i1 = 3;
P.i2 = 1;
% frequency grid
gridw = [0e0:2e-0:128];%129];%2.^[-2:0.1:10];%.1:1e-1:120;%

gridx = [0:1e0:50]*unit;
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
    ngy(:,t) = ngy(:,t)./sp(t);
    for i=1:5
        pfb(i,t) = sum(ngy(ifb{i},t));
    end
    mf(t) = sum(ngy(:,t).*gridw');
end
for i=1:5
    pfb(i,:) = pfb(i,:)./sum(pfb(i,:));
end

hf = figure('color',[1 1 1],'name','EEG frequency modulation');

ha = subplot(3,2,1,'parent',hf);
plot(ha,gridx,ds)
title(ha,'ds(z)/dz')
xlabel(ha,'fundamental mode dV (mV)')
set(ha,'xlim',[gridx(1),gridx(end)],'ygrid','on')
box(ha,'off')

ha = subplot(3,2,2,'parent',hf);
plot(ha,gridx,mf);
xlabel(ha,'fundamental mode dV (mV)')
title(ha,'centre frequency (Hz)')
set(ha,'xlim',[gridx(1),gridx(end)],'ygrid','on')
box(ha,'off')

ha = subplot(3,2,3,'parent',hf);
hm = mesh(ha,abs(gy).^2);
set(hm,'facecolor','flat','edgecolor',0.8*[1 1 1],'edgealpha',0.5)
axis(ha,'tight')
inf = find(gridw==round(gridw));
inf = inf(1:8:end);
set(ha,'ytick',inf,'yticklabel',gridw(inf),'xdir','reverse','zscale','log','cameraposition',[-88,472,0])
xlabel(ha,'fundamental mode dV (mV)')
ylabel(ha,'frequency (rad/s)')
title('frequency power')

ha = subplot(3,2,4,'parent',hf);
hm = mesh(ha,ngy);
set(hm,'facecolor','flat','edgecolor',0.8*[1 1 1],'edgealpha',0.5)
axis(ha,'tight')
% set(ha,'xticklabel',gridx(get(ha,'xtick')))
set(ha,'ytick',inf,'yticklabel',gridw(inf),'xdir','reverse','zscale','log','cameraposition',[-88,472,0])
xlabel(ha,'fundamental mode dV (mV)')
ylabel(ha,'frequency (rad/s)')
title('normalized frequency power')
colormap(flipud(autumn))

ha = subplot(3,2,5,'parent',hf);
plot(ha,gridx,sp);
xlabel(ha,'fundamental mode dV (mV)')
title(ha,'mean frequency power (A.U.)')
set(ha,'xlim',[gridx(1),gridx(end)],'ygrid','on')
box(ha,'off')

ha = subplot(3,2,6,'parent',hf);
plot(ha,gridx,pfb');
xlabel(ha,'fundamental mode dV (mV)')
title(ha,'band frequency power (A.U.)')
legend({'delta','theta','alpha','beta','gamma'})
set(ha,'xlim',[gridx(1),gridx(end)],'ygrid','on')
box(ha,'off')

getSubplots


f_fname = @f_modek;
g_fname = @g_Id;
n_t = 4e2;
u = zeros(1,n_t);
u(2) = 1;
options.inF = P;
options.inF.dt = 1e-3;
options.inF.k = [1:16];
options.inF.C = repmat([0;0;1],length(options.inF.k),1);
options.inF.z_10 = 0;
options.inG.ind = find(repmat([1;0;0],length(options.inF.k),1)==1);
options.dim = struct('n',3*length(options.inF.k),'n_phi',0,'n_theta',0,'n_t',n_t,'p',length(options.inF.k),'n_u',1);
x0 = repmat([0;0;0],length(options.inF.k),1);
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,[],[],u,Inf,Inf,options,x0);

displaySimulations(y,x,eta,e)





