function demo_WarOfAttrition

% demo for evolutionary 'war of attrition' game

% cf. 'Evolution and the theory of Games', Maynard-Smith, 1982.

V = 1; % value of the resource
k = 4; % scaling of quadratic time cost
quaCost = 1; % 0: linear cost in time, 1: quadratic cost in time
gridt = 0e-1:5e-2:4*V; % distribution of traits, in terms of costs

% Parameters of the simulation
n_t = 2e3;
dt = 1e-1;
f_fname = @log_replicator;
g_fname = @g_odds;
alpha   = Inf;
sigma   = Inf;
theta   = [V,k];
phi     = [];
u = [];


in.f_fitness = @fitness_woa;
in.dt = dt;
in.gridt = gridt;
in.quaCost = quaCost;
options.inF         = in;
options.inG         = in;
dim.n_theta         = length(theta);
dim.n_phi           = 0;
dim.n               = length(gridt);
dim.p               = length(gridt);
dim.n_t             = n_t;
options.dim = dim;

% Build time series of hidden states and observations
hf = figure('color',[1 1 1],'name','war of attrition');
ha = subplot(2,2,1,'parent',hf,'nextplot','add');
x0 = zeros(dim.n,1);
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);
imagesc(y,'parent',ha)
xlabel(ha,'time')
ylabel(ha,'cost prepared to pay')
title(ha,'phenotypes'' evolutionary dynamics')
axis(ha,'tight')
colorbar('peer',ha)
ha = subplot(2,2,2,'parent',hf);
imagesc(x,'parent',ha)
xlabel(ha,'time')
ylabel(ha,'cost prepared to pay')
title(ha,'deterministic log-odds dynamics')
axis(ha,'tight')
colorbar('peer',ha)

% evaluate stability of steady state
J = numericDiff(@f_replicator,1,y(:,end),theta,[],in);
ev = eig((J-eye(dim.n))./dt);
ha = subplot(2,2,3,'parent',hf,'nextplot','add');
rectangle('Position',[-1,-1,2,2],'curvature',[1 1],'parent',ha)
for i=1:length(ev)
    plot(ha,real(ev(i)),imag(ev(i)),'r+')
end
axis(ha,'equal')
grid(ha,'on')
xlabel(ha,'eigenvalues real axis')
ylabel(ha,'eigenvalues imaginary axis')
title(ha,'ESS: stability analysis')


% comparison of steady-state distribution and theoretical ESS
if ~quaCost
    p0 = (1/V).*exp(-gridt/V);
else
    p0 = (2*k*gridt/V).*exp(-k*gridt.^2/V);
end
p0 = p0./sum(p0);
p00 = y(:,end);
ha = subplot(2,2,4,'parent',hf,'nextplot','add');
plot(ha,gridt,p00,'r')
plot(ha,gridt,p0,'k--')
legend(ha,{'steady-state distribution','theoretical ESS'})
xlabel(ha,'cost prepared to pay')
ylabel(ha,'phenotypes'' frequency')
title(ha,'steady-state distribution versus theoretical ESS')

getSubplots

