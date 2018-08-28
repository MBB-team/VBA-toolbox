% Demo: on-line design optimization for DCM of fMRI data.
% This demo simplifies the network identification of Friston et al. (2003).
% In brief, photic input enters V1, which is reciprocally connected to V5.
% The goal of the experiment is to assess whether attention modulates the
% feedforward connection from V1 to V5. This is addressed using Bayesian
% model comparison given fMRI data (i.e. model 1: attention modulates the
% V1->V5 connection; model 2: no modulatory effect). Within the session,
% each block consists of 16 seconds of stimulation and 16 of rest. The
% on-line design optimization consists in deciding, before each block,
% whether we modulate subjects' attention. This is done by comparing the
% design efficiency  of two canonical block designs (i.e. with and without
% attentional modulation).
% NB: the procedure is adaptative, because it updates the posterior on
% unknown model parameters after each block. This is why the design
% efficiency inceases across time.



close all
clear variables


% Basic settings
n_t = 16; % number of time samples per block
nblocks = 10; % number of blocks
nreg = 2; % number of regions
TR = 3e0; % sampling period (in sec)
microDT = 5e-2; % micro-time resolution (in sec)
f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;


% construct atomic designs
nu = 2;
u1 = [0,ones(1,floor(n_t/2))];
u1(end:n_t) = 0;
u2 = zeros(1,n_t);
u = {[u1;u2];[u1;u1]}; % u{1}/u{2}: w/wo attentinal modulation

nm = 2; % # candidate models (see below)
truemodel = 1; % index of true model (which generates the data)


% model 1: u_att modulates connection from V1 to V5
A = [0 1
    1 0];
B{1} = zeros(nreg,nreg);
B{2} = zeros(nreg,nreg);
B{2}(2,1) = 1;
C = [1 0
     0 0];
D{1} = [0 0
    0 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);
[o4design{1},dim{1}] = getOptions4dcm(A,B,C,D,TR,microDT,n_t,1,1,1);


% model 1: no modulatory effect of attention
A = [0 1
    1 0];
nreg = size(A,1);
B{1} = zeros(nreg,nreg);
B{2} = zeros(nreg,nreg);
C = [1 0
     0 0];
D{1} = [0 0
    0 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);
[o4design{2},dim{2}] = getOptions4dcm(A,B,C,D,TR,microDT,n_t,1,1,1);


o4design{truemodel}.verbose = 0;

% simu parameters
t_A = exp([ -0.5
            -0.5
            ]);
t_Aself = -0;
t_B{1} = [];
t_B{2} = exp([ -0.5 ]);
t_C = exp([ -0.5 ]);
t_D{1} = [];
t_D{2} = [];
theta = zeros(dim{truemodel}.n_theta,1);
phi = zeros(dim{truemodel}.n_phi,1);
theta(o4design{truemodel}.inF.indA) = t_A;
for i=1:nu
    theta(o4design{truemodel}.inF.indB{i}) = t_B{i};
end
theta(o4design{truemodel}.inF.indC) = t_C;
alpha = Inf; % state noise precision
sigma = 5e-1; % measurement noise precision
x0 = zeros(dim{truemodel}.n,1); % initial conditions


% prepare graphics & montoring variables
hf = figure('color',[1 1 1]);
ha = subplot(3,2,1,'parent',hf,'nextplot','add','xlim',[1 nblocks+1]);
xlabel(ha,'time (blocks)')
ylabel(ha,'design efficiency')
ha2 = subplot(3,2,3,'parent',hf,'nextplot','add','xlim',[1 n_t*nblocks]);
xlabel(ha2,'time (scans)')
ylabel(ha2,'BOLD signal')
ha3 = subplot(3,2,2,'parent',hf,'nextplot','add','xlim',[1 nblocks+1]);
xlabel(ha3,'time (blocks)')
ylabel(ha3,'log p(y|m_1) -log p(y|m_2)')
e = zeros(2,1);
Y = zeros(dim{truemodel}.p,n_t*nblocks);
U = zeros(2,n_t*nblocks);
F = zeros(2,nblocks);
eb = zeros(1,nblocks);
vb = zeros(1,nblocks);

% on-line experiment

for tt=1:nblocks
    
    % 1- find best design given current information
    fprintf(1,'\n')
    fprintf(1,['-- Optimizing design (block ',num2str(tt),')...'])
    fprintf(1,'\n')
    for i=1:length(u) %loops over experimental designs
        fprintf(1,['- Candidate design #',num2str(i),'...'])
        fprintf(1,'\n')
        [e(i,tt)] = VBA_designEfficiency(f_fname,g_fname,dim,o4design,u{i},'models');
    end
    [em,ind] = max(e(:,tt));
    U(:,(tt-1)*n_t+1:tt*n_t) = u{ind};
    fprintf(1,['Optimizing design (block ',num2str(tt),')...  OK.'])
    fprintf(1,'\n')
    plot(ha,tt,e(1,tt)','ro')
    plot(ha,tt,e(2,tt)','go')
    legend(ha,{'u_{att} = off','u_{att} = on'},'Location','southeast','Orientation','horizontal')
    
    % 2- simulate BOLD response to chosen design under true model
    fprintf(1,['-- Simulating data (block ',num2str(tt),')...  '])
    [y,x,x0,eta,ee] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u{ind},alpha,sigma,o4design{truemodel},x0);
    Y(:,(tt-1)*n_t+1:tt*n_t) = y;
    x0 = x(:,end);
           
    try
        set(pl(1), 'XData', 1 : tt*n_t);
        set(pl(1), 'YData', [get(pl(1), 'YData') y(1,:)]);
        set(pl(2), 'XData', 1 : tt*n_t);
        set(pl(2), 'YData', [get(pl(2), 'YData') y(2,:)]);
    catch
        pl(1) = plot(ha2,1:n_t,y(1,:),'m');
        pl(2) = plot(ha2,1:n_t,y(2,:),'b');
    end  
    legend(ha2,{'V1','V5'},'Location','southeast','Orientation','horizontal')
    fprintf(1,[' OK.'])        
    fprintf(1,'\n')
    drawnow
    
    % 3- invert both models given new piece of dataset
    % NB: priors for the next block are updated to current posterior
    fprintf(1,['-- VB (block ',num2str(tt),'):   inverting model       '])
    for j=1:length(o4design)
        fprintf(1,repmat('\b',1,6))
        fprintf(1,[num2str(j),'/',num2str(length(o4design)),'...'])
        o4design{j}.DisplayWin = 0;
        o4design{j}.verbose = 0;
        [posterior,out] = VBA_NLStateSpaceModel(y,u{ind},f_fname,g_fname,dim{j},o4design{j});
        o4design{j}.priors = posterior;
        o4design{j}.priors.muX0 = posterior.muX(:,end);
        o4design{j}.priors.SigmaX0 = posterior.SigmaX.current{end};
        if tt > 1
            F(j,tt) = out.F + F(j,tt-1);
        else
            F(j,tt) = out.F;
        end
        OUT(j,tt).out = out;
        OUT(j,tt).posterior = posterior;
    end
    plot(ha3,tt,F(1,tt)-F(2,tt),'k*')
    
    eb(tt) = OUT(1,tt).posterior.muTheta(OUT(1,tt).out.options.inF.indB{2});
    vb(tt) = OUT(1,tt).posterior.SigmaTheta(OUT(1,tt).out.options.inF.indB{2},OUT(1,tt).out.options.inF.indB{2});
    
    fprintf(1,repmat('\b',1,24))
    fprintf(1,[' OK.'])        
    fprintf(1,'\n')
    
    
end

ha4 = subplot(3,2,4,'parent',hf,'nextplot','add','xlim',[0 nblocks+1]);
xlabel(ha4,'time (blocks)')
ylabel(ha4,'modulatory effect')
plotUncertainTimeSeries(eb',1.96^2*vb',[],ha4)
hold(ha4,'on')
plot(ha4,[0,nblocks+1],[t_B{2},t_B{2}],'g--')

ha5 = subplot(3,2,5,'parent',hf);
imagesc(U,'parent',ha5)
title(ha5,'chosen (online) design')
xlabel(ha5,'time (scans)')
hold(ha5,'on')
plot(get(ha5,'xlim'),[1.5 1.5],'k')
set(ha5,'ytick',[1,2],'yticklabel',{'u1','u2'})
colormap(flipud(bone))

% invert full datasets at once
n_t = size(Y,2);
B{2}(2,1) = 1; % add modulatory effect (true model)
[OPT,DIM] = getOptions4dcm(A,B,C,D,TR,microDT,n_t,1,1,1);
[posterior,out] = VBA_NLStateSpaceModel(Y,U,f_fname,g_fname,DIM,OPT);

