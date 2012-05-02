% Demo: design optimization for DCM of fMRI data.


close all
clear variables



%-----------------------------------------------------------
%-------------- DCM model specification --------------------

%--------- Basic settings ----------------------------------
n_t = 11;                    % number of time samples
TR = 3e0;                       % sampling period (in sec)
microDT = 5e-2;                 % micro-time resolution (in sec)
f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;





%--------- model 1 -----------------------------------------
% invariant effective connectivity
A = [0 1
    1 0];
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
% input-state coupling
C = [1;0];
% gating (nonlinear) effects
D{1} = [0 0
    0 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);
% Build priors and options/dim structures for model inversion
[o4design{1},dim{1}] = getOptions4dcm(A,B,C,D,TR,microDT,n_t);


%--------- model 2 -----------------------------------------
% invariant effective connectivity
A = [0 1
    1 0];
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
% input-state coupling
C = [0;1];
% gating (nonlinear) effects
D{1} = [0 0
    0 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);
% Build priors and options/dim structures for model inversion
[o4design{2},dim{2}] = getOptions4dcm(A,B,C,D,TR,microDT,n_t);


u0 = zeros(1,n_t);

nm = length(dim);

truemodel = 1;


%---------- Simu parameters
t_A = exp([ -0.5
            -0.5
            ]);
t_Aself = -0;
t_B{1} = [];
t_C = exp([ +0.1 ]);
t_D{1} = [];
t_D{2} = [];
theta = zeros(dim{truemodel}.n_theta,1);
phi = zeros(dim{truemodel}.n_phi,1);
theta(o4design{truemodel}.inF.indA) = t_A;
theta(o4design{truemodel}.inF.indC) = t_C;
alpha = Inf;               % state noise precision
sigma = 1e1;                  % measurement noise precision




nblocks = 10;

hf = figure;
ha = axes('parent',hf,'nextplot','add');

e = zeros(2,nblocks);

for tt=1:nblocks
    
    
    u{1} = u0;
    u{1}(end-10:end-5) = 1;
    
    u{2} = u0;
    
    
    for i=1:length(u) %loops over experimental designs
        for j=1:length(dim)
            dim{j}.n_t = length(u{i}(end-10:end));
        end
        
        [e(i,tt)] = designEfficiency(f_fname,g_fname,dim,o4design,u{i}(end-10:end),'models');
        

        
    end
    
    plot(ha,[tt;tt],e(:,tt),'o')
%     pause
    
    [em,ind] = max(e(:,tt));
    uu = u{ind};
    
    [y,x,x0,eta,ee] = simulateNLSS(length(uu),f_fname,g_fname,theta,phi,uu,alpha,sigma,o4design{truemodel});
    displaySimulations(y,x,eta,ee)
%     pause
    
    
    for j=1:length(o4design)
        dim{j} = rmfield(dim{j},'n_t');
        o4design{j}.DisplayWin = 0;
        o4design{j}.verbose = 0;
        [posterior,out] = VBA_NLStateSpaceModel(y(:,end-10:end),uu(:,end-10:end),f_fname,g_fname,dim{j},o4design{j});
        o4design{j}.priors = posterior;
        o4design{j}.priors.muX0 = posterior.muX(:,end);
        o4design{j}.priors.SigmaX0 = posterior.SigmaX.current{end};
    end
    
%     pause
    
    u0 = [uu,zeros(1,n_t)];
    
    
    
    
end


