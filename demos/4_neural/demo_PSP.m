% demo for postsynaptic potentials simulations

function demo_PSP

n_t = 1e3;
u = zeros(1,n_t);
u(2) = 1;

theta = [1;5e-3;150e-3];
% theta(1): amplitude of the alpha kernel
% theta(2): rising time
% tehta(3): decay time
phi = [];

f_fname = @f_PSP;
g_fname = @g_Id;

dim.n = 2;
dim.n_theta = 2;
dim.n_phi = 0;

options.inG.ind = 1;
options.inF.dt = 1e-3;

x0 = [0;0];

[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0);

displaySimulations(y,x,eta,e);



