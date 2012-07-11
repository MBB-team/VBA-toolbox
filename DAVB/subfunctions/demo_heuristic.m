% demo for the Kilner heuristic

% copute dispersive propagator
gt = 1e-2:1e-4:2e-1;
gr = -4e-2:1e-4:4e-2;
n1 = length(gt);
n2 = length(gr);
G = zeros(n1,n2);
P = [3e-1;3e-3;2e3];
for i=1:n2
    [G(:,i)] = dispersivePropagator(gr(i),gt(:),P);
end
figure,imagesc(G)
[X,Y] = meshgrid(gt,gr);
figure,
h = mesh(X,Y,G');
i0 = find(gr==0);
i1 = find(gt==2e-2);
hold on
plot3(gt,0.*gt,G(:,i0)','color',[1 1 1])
plot3(gt(i1)*ones(size(gr)),gr,G(i1,:)','color',[1 1 1])
plot3(gt,5e-2*ones(size(gt)),G(:,i0)','k--')
hp = plot3(-2e-2.*ones(size(gr)),gr,G(1,:),'k--')
set(gca,'ytick',[-4e-2:2e-2:4e-2],'yticklabel',[-4:2:4])
set(gca,'xtick',[0:5e-2:2e-1],'xticklabel',[0:50:200])
set(gca,'ylim',[-4e-2,5e-2],'xlim',[-2e-2,2e-1])
xlabel('time (msec)')
ylabel('distance (mm)')

return

% excitatory maximum post-synaptic depolarization
P.me = 1e-0*8;
% inhibitory maximum post-synaptic depolarization
P.mi = 1e-0*32;
% excitatory post-synaptic time constant
P.Ke = 1e-3/4;
% inhibitory post-synaptic time constant
P.Ki = 1e-3/28;
% amplitude of intrinsic connectivity kernels
P.a = 1e3.* [ 0 0 2
              0 0 8
              2 1 0 ];
% intrinsic connectivity decay constant
P.c = 1e3*0.32.*ones(3,3); 
% conduction velocity
P.v = 3.*ones(3,3); 
% radius of cortical source
P.l = 50.*1e-3;
% parameters of the Gaussian observation filter
P.phi = [1;0.01.*P.l];
% parameters of the Gaussian observation filter
P.sig = struct('r',0.54,'eta',30*1e-0,'g',0.135);
% frequency grid
gridw = .1:1e-1:120;

gridx = -50:1e0:50;
nt = numel(gridx);
gy = zeros(numel(gridw),nt);

for t=1:nt
    [g] = dsdv(gridx(t),P);
    [gy(:,t)] = spectralPower2(P,gridw,g)';
    t
end

ngy = gy;
mf = zeros(1,nt);
for t=1:nt
    ngy(:,t) = ngy(:,t)./sum(ngy(:,t));
    mf(t) = sum(sqrt(ngy(:,t)).*gridw');
end

hf = figure;
ha = subplot(2,1,2,'parent',hf);
plot(ha,gridx,mf);
ha = subplot(2,1,1,'parent',hf);
mesh(ha,log(gy));
axis(ha,'tight')


