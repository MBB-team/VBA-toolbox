function [posterior_fullW,out_fullW,ctrst] = demo_negfeedback(As, Bs, hAs, noise,reps)
%%% ----------------------------------------------------------------------
%   Initialise the generative model, simulate data and inverse the model
%   
%   IN:     feedback 
%%% ----------------------------------------------------------------------

% close all

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

if nargin==0
   As = [1 0];
   Bs = [0 -1.5] ;
   hAs = [1 0];  
   noise=1;
   reps=3;
end

% === Basic settings =======================================
f_fname = @f_DCMwHRFext ;  %
g_fname = @g_demo_extended;   % 

                              % 
TR = .5;                      % sampling period (in sec)
n_t = round(160/TR);          % number of time samples
microDT = .05;               % micro-time resolution (in sec)
homogeneous = 0;              % params of g(x) homogeneous accross regions
reduced_f = 0;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;                  % state noise precision
sigma = 1/noise;                % measurement noise precision

% === Input ================================================

u=[];
answers = [];
ctrst=[];

isi = zeros(2,15) ;
u0 = 0*ones(2,15);
ctrst0 = 0*ones(2,15);
answers0 = zeros(1,15);

for u1=[1]
    for u2=[1 .01]
        u_temp = u0;
        u_temp(1,1:2) = u1 ;
        u_temp(2,:) = u2 ;
        ctrst0(1,2:6) = u1;
        ctrst0(2,2:6) = u2;

        answers_temp=answers0;
        answers_temp(2:6)=1:5;
        
%         jitter = zeros(2,randi(6)-1);
        u=[u isi  u_temp]; 
        ctrst=[ctrst isi  ctrst0];
        answers=[answers isi(1,:)  answers_temp];
   
    end
end

% reps=5;
u = [repmat(u,1,reps) isi];
answers = [repmat(answers,1,reps) isi(1,:)];
ctrst = [repmat(ctrst,1,reps) isi];


nu = size(u,1);
n_t=size(u,2);


% === DCM structure ========================================
% invariant effective connectivity
A = [0 1;
     1 1];

nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = [0 1;
        1 0];
% input-state coupling
C = [1 0;
     0 0];
% gating (nonlinear) effects
D{1} = zeros(nreg,nreg);
D{2} = zeros(nreg,nreg);
% === Decoding scheme ========================================
hA = [1 1];
nrep=size(hA,1);
hB = {zeros(nrep,nreg),zeros(nrep,nreg)}; 
hC = zeros(nrep,nu);
hD = {  zeros(nrep,nreg),...
        zeros(nrep,nreg) ...
     };

% === Build options and dim structures =====================
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);

options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = 8;
options.GnFigs = 0;
options.inF.linearized = lin;
options.DisplayWin=1;
options.verbose=1;
dim.n_theta = options.inF.ind5(end);
if options.extended
    dim.n_phi = options.inG.indr;
else
    dim.n_phi = options.inG.ind2(end);
end
dim.n = 5*nreg+nrep;
dim.n_t=n_t;


sources(1) = struct('out',1:2,'type',0); % BOLD signal (gaussian, dim=4)
sources(2) = struct('out',3,'type',1); % first binary response
options.sources=sources;


%% -----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
theta = zeros(dim.n_theta,1);
phi = zeros(dim.n_phi,1);

%- DCM
t_Aself = log(1);


t_A = [As .4]; %[1.5 -1 .15];

t_B{1} = [];
t_B{2} = Bs ; 
t_C = [1];
t_D{1} = [];
t_D{2} = [];

theta(options.inF.indA) = t_A;
theta(options.inF.indself) = t_Aself;
for i=1:nu
    theta(options.inF.indB{i}) = t_B{i};
end
theta(options.inF.indC) = t_C;
for i=1:nreg
    theta(options.inF.indD{i}) = t_D{i};
end
%- decoding 

t_hA = hAs*1.5; 

t_hB{1} = [];
t_hB{2} = [];
t_hC = [];
t_hD{1} = [];
t_hD{2} = [];
t_hself = log(1/options.inF.deltat) ;

theta(options.inF.indhA) = t_hA;
theta(options.inF.indhself) = t_hself;
for i=1:nu
    theta(options.inF.indhB{i}) = t_hB{i};
end
theta(options.inF.indhC) = t_hC;
for i=1:nreg
    theta(options.inF.indhD{i}) = t_hD{i};
end

%- observation
phi(options.inG.indr) = 0;

%--- Simulate time series of hidden states and observations
%disp('*** Simulation');
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
% plot(y(1:2,:)');
% hold on;
% temp=y(3,:);
% temp(answers==0)=NaN;
% plot(temp,'rx');
% hold off;

perfs=[];
us=unique(u(2,u(2,:)>0));
    for ui=1:length(us)
        for rt=unique(answers(answers>0))
            perfs(ui,rt) = mean(y(3,find(answers==rt & u(2,:)==us(ui))));
        end
        if us(ui)==1
            perfs(ui,:) = perfs(ui,:);
            perfs(ui+1,:) = 1-perfs(ui,:);
        end
    end

% figure();
% plot(perfs');

% display time series of hidden states and observations
% displaySimulations(y,x,eta,e)
%  disp('--paused--')
%  pause
% figure()
% subplot(2,1,1)
% plot(perfs');hold on;plot(1:size(perfs,2),.5*ones(1,size(perfs,2)),'k:'); hold off
% 
% subplot(2,1,2)
% temp=y(3,1:180)-e(3,1:180);plot(temp);temp(answers(1:length(temp))==0)=NaN;hold on; plot(temp,'.');ylim([0,1]);hold off;
% 
% with=[];
% without=[];
% return
%-----------------------------------------------------------
%------------------- model inversion -----------------------
%--- Call inversion routine
% disp('*** Hypothesis inversion');

% -- decoded model
hA = [1 1];
A = [0 1;
     1 1];
 B{2} = [0 1;
        1 0];
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);

options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = 8;
options.GnFigs = 0;
options.inF.linearized = lin;
options.DisplayWin=1;
options.verbose=1;
options.graphic=1;
options.checkGrads = 0;
dim.n_theta = options.inF.ind5(end);
if options.extended
    dim.n_phi = options.inG.indr;
else
    dim.n_phi = options.inG.ind2(end);
end
dim.n = 5*nreg+nrep;
dim.n_t=n_t;


sources(1) = struct('out',1:2,'type',0); % BOLD signal (gaussian, dim=4)
sources(2) = struct('out',3,'type',1); % first binary response
options.sources=sources;

options.isYout=zeros(size(y));
options.isYout(1,:)=mod(1:size(y,2),4)>0;
options.isYout(2,:)=options.isYout(1,:);
options.isYout(3,:)=1*(answers==0);

[posterior_fullW,out_fullW] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior_fullW,out_fullW,y-e,x,x0,theta,phi,alpha,sigma)

return

prior_with = options.priors;
dim_with=out_fullW.dim;

% -- fMRI only model model
options.isYout=zeros(size(y));
options.isYout(1,:)=mod(1:size(y,2),4)>0;
options.isYout(2,:)=options.isYout(1,:);
options.isYout(3,:)=1;
% no need of decoding scheme
options.priors.SigmaTheta(options.inF.indhA,options.inF.indhA)=0;

[posterior_fullWO,out_fullWO] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

prior_without = options.priors;
dim_without=out_fullWO.dim;

Fs_with=[];
Fs_without=[];
dscheme = {[1 0]; [0 1]} ;%; 1 1];
nested= {[1 0 0 1]; [0 1 1 0] ; [1 1 0 0]};
id_netsted = [options.inF.indA(1:2) options.inF.indB{2}];
id_scheme = options.inF.indhA ;
for nested_i = 1:length(nested) 
    for dscheme_i = 1:2
    priors2_with = prior_with;
    priors2_without = prior_without;
    
    priors2_with.SigmaTheta(id_netsted,id_netsted)=diag(nested{nested_i}).*prior_with.SigmaTheta(id_netsted,id_netsted);
    priors2_without.SigmaTheta(id_netsted,id_netsted)=diag(nested{nested_i}).*prior_without.SigmaTheta(id_netsted,id_netsted);
   
    priors2_with.SigmaTheta(id_scheme,id_scheme)=diag(dscheme{dscheme_i})*prior_with.SigmaTheta(id_scheme,id_scheme);
    
    [F_nestedW,~] = VB_SavageDickey(posterior_fullW,prior_with,out_fullW.F,dim_with,priors2_with) ;
	[F_nestedWO,~] = VB_SavageDickey(posterior_fullWO,prior_without,out_fullWO.F,dim_without,priors2_without) ;
    
    Fs_with(dscheme_i,nested_i) = F_nestedW;
    Fs_without(dscheme_i,nested_i) = F_nestedWO;
    end
    
end


with.connections = posterior_fullW.muTheta;
with.F = out_fullW.F;
with.Fnested =Fs_with;

without.connections = posterior_fullWO.muTheta;
without.F = out_fullWO.F;
without.Fnested =Fs_without;



data.posterior_fullW=posterior_fullW;
data.out_fullW=out_fullW;

end