% demo for dynamical systems simulations

close all
clear variables

% % simulate simple 1D linear system
% f_fname = @f_L1;
% g_fname = @g_Id;
% n_t = 1e2;            % number of time samples
% dt = 2e-2;               % micro-time resolution (in sec)
% alpha = Inf;              % state noise precision
% sigma = Inf;              % measurement noise precision
% u = zeros(1,n_t);
% % u(1) = 1;
% dim.n_theta = 2;
% dim.n_phi = 0;
% dim.n = 1;
% options.inF.dt = dt;
% N = 16;
% hf = figure('color',[1 1 1]);
% ha = subplot(1,2,1,'parent',hf,'nextplot','add');
% ha(2) = subplot(1,2,2,'parent',hf,'nextplot','add');
% for i=1:N
%     x0 = 4*randn;
%     theta = [-1.2;1];
%     [y,x] = simulateNLSS(n_t,f_fname,g_fname,theta,[],u,alpha,sigma,options,x0);
%     plot(ha(1),[0:n_t-1].*dt,x)
%     plot(ha(1),0,x0,'.')
%     set(ha(1),'box','off','ygrid','on','xlim',[-0.5;2],'ylim',[-22;22]);
%     xlabel(ha(1),'time (sec)')
%     ylabel(ha(1),'x(t)')
%     title(ha(1),['a=',num2str(theta(1),2)])
%     theta = [1.2;1];
%     [y,x] = simulateNLSS(n_t,f_fname,g_fname,theta,[],u,alpha,sigma,options,x0);
%     plot(ha(2),[0:n_t-1].*dt,x)
%     plot(ha(2),0,x0,'.')
%     set(ha(2),'box','off','ygrid','on','xlim',[-0.5;2],'ylim',[-22;22]);
%     xlabel(ha(2),'time (sec)')
%     ylabel(ha(2),'x(t)')
%     title(ha(2),['a=',num2str(theta(1),2)])
% end
% plot(ha(1),get(ha(1),'xlim'),[0,0],'r--')
% plot(ha(2),get(ha(2),'xlim'),[0,0],'r--')
% getSubplots



% simulate alpha-kernel
f_fname = @f_alpha;
g_fname = @g_Id;
n_t = 1e3;            % number of time samples
dt = 1e-1;               % micro-time resolution (in sec)
alpha = Inf;              % state noise precision
sigma = Inf;              % measurement noise precision
u = zeros(1,n_t);
times = ceil(n_t*rand(16,1));
u(times) = 1;
dim.n_theta = 2;
dim.n_phi = 0;
dim.n = 1;
options.inF.dt = dt;

theta = [1/dt;1];
x0 = [0;0];
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,[],u,alpha,sigma,options,x0);

hf = figure('color',[1 1 1]);
ha = subplot(1,2,1,'parent',hf,'nextplot','add');
plot(ha(1),[0:n_t-1].*dt,u)
set(ha(1),'box','off','ygrid','on','ylim',[-0.2,1.2]);
xlabel(ha(1),'time (sec)')
ylabel(ha(1),'u(t)')
title(ha(1),['input u'])
ha(2) = subplot(1,2,2,'parent',hf,'nextplot','add');
plot(ha(2),[0:n_t-1].*dt,x(1,:))
plot(ha(2),0,x0,'.')
set(ha(2),'box','off','ygrid','on','ylim',[-0.2,1.2]);
xlabel(ha(2),'time (sec)')
ylabel(ha(2),'x(t)')
title(ha(2),['output x'])

% show likelihood
yg = [-0.2:1e-2:1.2];
LL = zeros(length(yg),n_t);
for t=1:n_t
    LL(:,t) = exp(-32*(yg(:)-y(1,t)).^2);
end
hf = figure('color',[1 1 1]);
ha = subplot(1,2,1,'parent',hf,'nextplot','add');
imagesc(LL,'parent',ha)
set(ha(1),'box','off','xlim',[0,n_t-1]);%,'ylim',[-0.2,1.2]);
axis(ha(1),'off')
xlabel(ha(1),'time')
ylabel(ha(1),'p(y|P,m)')
title(ha(1),['likelihood'])

gp1 = theta(1)+[-4:4];
gp2 = theta(2)+[-2:0.5:2];
n = length(gp);
for i=1:n
    i
    for j=1:n
        theta = [gp1(i);gp2(j)];
        [gy] = simulateNLSS(n_t,f_fname,g_fname,theta,[],u,alpha,sigma,options,x0);
        LLp(i,j) = sum((gy(1,:)-y(1,:)).^2);
    end
end
LLp(find(isnan(LLp)==1) = Inf;
LLp(LLp>1e3) = 1e3;
ha(2) = subplot(1,2,2,'parent',hf,'nextplot','add');
imagesc(exp(-1e-2*LLp),'parent',ha(2))
set(ha(2),'box','off');
axis(ha(2),'off')
xlabel(ha(2),'time')
ylabel(ha(2),'p(y|P,m)')
title(ha(2),['likelihood'])
      

options.priors.muTheta = .5*ones(2,1);
dim.n = 2;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);



