% Demo for DCM for fMRI: motor and premotor interactions
% This demo looks for motor/premotor responses to 'prepare' and 'go'
% instructions.

close all
clear variables

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

%--- Basic settings
f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;
TR = 2e0;                     % sampling period (in sec)
T = 2e2;                     % total sampling time (in sec)
n_t = round(T/TR);          % number of time samples             
microDT = 2e-1;               % micro-time resolution (in sec)
homogeneous = 0;              % params of g(x) homogeneous accross regions
reduced_f = 0;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;%1e2/TR;               % state noise precision
sigma = 1e0;                  % measurement noise precision

%--- Input
iu_prepare = 10+[0:50:150]./microDT;
iu_go = iu_prepare(1:2:end)+10/microDT;
u_prepare = zeros(1,T/microDT);
u_go = zeros(1,T/microDT);
for i=1:length(iu_prepare)
    u_prepare(iu_prepare(i):iu_prepare(i)+10/microDT-1) = 1;
end
for i=1:length(iu_go)
    u_go(iu_go(i):iu_go(i)+2) = 1;
end
u(1,:) = u_prepare+u_go;
u(2,:) = u_go;

nu = size(u,1);

%--- DCM structure: model 'drive'
% invariant effective connectivity
A = [0 1
     1 0];
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = [0 0
        1 0];
% input-state coupling
C = [1 0
     0 0];
% gating (nonlinear) effects
D{1} = zeros(nreg,nreg);
D{2} = zeros(nreg,nreg);


%--- Build options and dim structures
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous);
options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 1;
options.backwardLag = ceil(16/TR);  % 16 secs effective backward lag
options.inF.linearized = lin;
dim.n_theta = options.inF.ind5(end);
if options.inG.homogeneous
    dim.n_phi = 2;
else
    dim.n_phi = 2*nreg;
end  
dim.n = 5*nreg;


%-----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
t_A = [0;0];
t_Aself = -0;
t_B{1} = [];
t_B{2} = 4;
t_C = exp([ +0.1 ]);
t_D{1} = [];
t_D{2} = [];
theta = zeros(dim.n_theta,1);
phi = zeros(dim.n_phi,1);
theta(options.inF.indA) = t_A;
for i=1:nu
    theta(options.inF.indB{i}) = t_B{i};
end
theta(options.inF.indC) = t_C;



%--- Simulate time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
displaySimulations(y,x,eta,e);
figure,imagesc(u)
figure,plot(x(options.inF.n5,:)')
% pause


% %--- Call inversion routine
% [p1,o1] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
% 
% %--- Display results
% displayResults(p1,o1,y-e,x,x0,theta,phi,alpha,sigma)
% 
% 






%--- DCM structure: model 'withold'

u(1,:) = u_prepare;
u(2,:) = u_prepare+u_go;
u(3,:) = 1;

% invariant effective connectivity
A = [0 1
     1 0];
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = zeros(nreg,nreg);
B{3} = [1 0
        0 1];
% input-state coupling
C = [1 0 0
     0 1 0];
% gating (nonlinear) effects
D{1} = zeros(nreg,nreg);
D{2} = zeros(nreg,nreg);


%--- Build options and dim structures
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous);
options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 1;
options.backwardLag = ceil(16/TR);  % 16 secs effective backward lag
options.inF.linearized = lin;
dim.n_theta = options.inF.ind5(end);
if options.inG.homogeneous
    dim.n_phi = 2;
else
    dim.n_phi = 2*nreg;
end  
dim.n = 5*nreg;

options.priors.muTheta(options.inF.indB{3}(1)) = -5;

%--- Call inversion routine
[p2,o2] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
set(gcf,'name','model ''withold''')
VBA_spm_dcm_explore(vba2dcm(p2,o2,[],TR));


% VBA_ReDisplay(p1,o1,1)
% set(gcf,'name','model ''drive''')
