function [posterior,out,simulated] = demo_behaviouralDCM()
%%% ----------------------------------------------------------------------
%%% ----------------------------------------------------------------------

% close all

%-----------------------------------------------------------
%-------------- DCM model specification --------------------


% === Basic settings =======================================
f_fname = @f_DCMwHRFext ;  %
g_fname = @g_DCMwHRFext;   % 

                              % 
TR = 2;                      % sampling period (in sec)
microDT = .1;               % micro-time resolution (in sec)

alpha = Inf;                  % state noise precision
sigma = 5e1;                % measurement noise precision

% === Input ================================================

repets = 20;
u = repmat( ...
    [0  1 0 0 0 0 0  1 0 0 0 0 ;
     0 +1 0 0 0 0 0 -1 0 0 0 0] ...
    ,1,repets);

[nu,n_t]=size(u);

isResponse = u(1,:) == 1;

% === DCM structure ========================================
% invariant effective connectivity
A = [0 0;
     0 0];

nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = [0 0;
        1 0];
% input-state coupling
C = [1 0;
     0 0];
% gating (nonlinear) effects
D{1} = zeros(nreg,nreg);
D{2} = zeros(nreg,nreg);

% === Decoding scheme ========================================
hA = [0 1];
nrep=size(hA,1);
hB = {zeros(nrep,nreg),zeros(nrep,nreg)}; 
hC = zeros(nrep,nu);
hD = {  zeros(nrep,nreg),...
        zeros(nrep,nreg) ...
     };

% === Build options and dim structures =====================
sources(1) = struct('out',1:2,'type',0); % BOLD signal (gaussian, dim=4)
sources(2) = struct('out',3,'type',1); % first binary response

options = prepare_fullDCM(A,B,C,D,TR,microDT,1,hA,hB,hC,hD,sources);
options.dim.n_t = n_t;
dim = options.dim ;


%% -----------------------------------------------------------
%----------- simulated times series specification ----------

theta = zeros(dim.n_theta,1);
phi   = zeros(dim.n_phi  ,1);

%- DCM
t_Aself = log(10);
t_A = []; 

t_B{1} = [];
t_B{2} = [3] ;

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
t_hA = [2]; 

t_hB{1} = [];
t_hB{2} = [];
t_hC = [];
t_hD{1} = [];
t_hD{2} = [];

t_hself = log(1) ;

theta(options.inF.indhA) = t_hA;
theta(options.inF.indhself) = t_hself;
for i=1:nu
    theta(options.inF.indhB{i}) = t_hB{i};
end
theta(options.inF.indhC) = t_hC;
for i=1:nreg
    theta(options.inF.indhD{i}) = t_hD{i};
end

simulated.theta = theta;
simulated.phi = phi;

%--- Simulate time series of hidden states and observations
disp('*** Simulation');

[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,zeros(dim.n,1));


f=figure('Color','w');
subplot(2,1,1)
plot((y(1:2,:)-e(1:2,:))'); hold on
plot(y(1:2,:)','.');                
ylim([-.8 1.6]);
ylabel('BOLD signal')

subplot(2,1,2)
plot((y(3,:)-e(3,:))','r'); hold on
plot(find(isResponse),y(3,find(isResponse))','r.');          
ylim([-.1 1.1]);
ylabel('response')


options.isYout(3,:)=~isResponse;

%-----------------------------------------------------------
%------------------- model inversion -----------------------
disp('*** Inversion');

options.priors = getPriors(nreg,n_t,options,1,0);
[options.priors.a_sigma, options.priors.b_sigma] = VBA_guessHyperpriors(y(1:2,:),[0.05, 0.25]) ;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

set(0,'CurrentFigure',f);
subplot(2,1,1)
plot(out.suffStat.gx(1:2,:)',':'); hold off

subplot(2,1,2)
plot(out.suffStat.gx(3,:)',':r'); hold off


