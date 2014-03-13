function demo_evol_lossAversion

% demo for evolutionary selection of loss aversion

% Loss-aversion trait is defined on a 1D continuous interval, which is the
% range of a parameter controlling the gain/loss asymmetry of the utility
% function of agents (see function lossAverseUtile.m). When this parameter
% is equal to zero (resp., one), agents are insensitive to gains (resp.,
% losses).
% Let lambda be the parameter controlling the gain/loss asymmetry, and
% x(lambda) be the frequency of this trait within the population. This demo
% uses standard replicator dynamics of evolutionary game theory to derive
% the adaptive fitness of the trait, ie:
%
% dx_i/dt = x_i ( f_i(x) + sum_j{x_j*f_j(x)} )
%
% where f_i(x) is the malthusian fitness of trait i (ie of a given lambda).
% This fitness is obtained according to the expected outcome of a game
% (here: "stag-hunt"), for each type of agent or trait.
% In our case, lambda changes the probability to cooperate, through a
% probabilistic softmax action emission law.
%
% The outcome table of the "stag-hunt" game is as follows:
% A = [ X(1)-c  X(3)-c
%       X(2)-c  X(2)-c  ]
% where the first line of A applies when the player cooperates, and the
% second line applies when she defects, rather than colmuns pertain to the
% opponent's behaviour. Note that the payoff of defecting is independent of
% the opponent's behaviour.
% There is an effort cost (c) attached to any action. This cost does not
% affect the action emission law, since the difference in values is
% invariant to c. However, it turns out to be critical for deriving the
% evolutionary fitness of loss aversion.
% Note: the malthusian fitness of player i is its payoff, averaged across
% the population. The chance of encountering any member of type j is
% proportional to its frequency wihin the population. This means the rate
% of change of trait's frequencies within the population is a function of
% both the expected payoff (given expected actions of each type of agents)
% and the current distribution of trait frequencies.

n = 32; % number of bins for gridding the asymmetry parameter lambda.

% Parameters of the simulation
n_t = 2e3;
dt = 1e0;
f_fname = @log_replicator;
g_fname = @g_odds;
alpha   = Inf;
sigma   = Inf;
theta   = [];
phi     = [];
u = [];

in.X = [3;1;0]*1e0;
in.mc = 2; % mean cost (of effort)
in.sc = 0;
in.dt = dt;
in.lambda = 0:1./(n-1):1;
in.f_fitness = @fitness_lossAversion;


options.inF         = in;
options.inG         = in;
dim.n_theta         = 0;
dim.n_phi           = 0;
dim.n               = n;

% Build time series of hidden states and observations
x0 = zeros(n,1);

options.checkGrads = 0;

[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

% evaluate stability of steady state
J = numericDiff(@f_replicator,1,y(:,end),[],[],in);
ev = eig((J-eye(dim.n))./dt);
hf = figure('color',[1 1 1],'name','ESS: stability analysis');
ha = axes('parent',hf,'nextplot','add');
rectangle('Position',[-1,-1,2,2],'curvature',[1 1],'parent',ha)
for i=1:length(ev)
    plot(ha,real(ev(i)),imag(ev(i)),'r+')
end
axis(ha,'equal')
grid(ha,'on')

% display results
[my,mi] = max(y(:,end));
xg = -10:0.01:10;
Umap = lossAverseUtile(xg,in.lambda(mi));
Umean = lossAverseUtile(xg,sum(y(:,end).*in.lambda'));
U0 = lossAverseUtile(xg,1/2);
figure,plot(xg,Umap)
hold on
plot(xg,Umean,'r')
plot(xg,U0,'k--')
legend({'lambda: map','lambda: mean','lambda:neutral'})
xlabel('outcome numerical value')
ylabel('outcome utility')

figure,imagesc(y)
[hi,px,gx,x] = plotGraph3D(y',in.lambda');



gridc = [-4:0.2:4];
nc = length(gridc);
Y = zeros(size(y,1),nc);
for i=1:nc
    in.mc = gridc(i);
    options.inF         = in;
    options.inG         = in;
    [y] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);
    Y(:,i) = y(:,end);
end
figure,imagesc(Y)
[hi,px,gx,x] = plotGraph3D(Y,gridc');
Y;



