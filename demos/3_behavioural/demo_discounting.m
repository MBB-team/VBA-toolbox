% This demo compares exponential vs hyperbolic temporal discounting.
% This script simulates an agent making choices between two options, which
% differ in terms of the expected reward and temporal horizons, under both
% an hyperbolic and an exponantial temporal discounting model.
% The choice data are then inverted using these two models. NB: in both
% models, the discount factor and the behavioural temperature are unknown.
% The effect of simulating the data under each model onto the log-Bayes
% factor (comparison between hyperbolic vs exponential models) is assessed
% using Monte-Carlo simulations. 

clear all
close all

models = {'hyperbolic','exponential'};

% initialize simulations
ntrials = 128;
g_fname = @g_discounting;
phi = [1;-1];
in.ind.logk = 1; % index of discount factor (phi)
in.ind.logb = 2; % index of behavioural temperature (phi)
in.ind.t = [1;3]; % index of temporal horizon (u)
in.ind.R = [2;4]; % index of expected reward (u)



% Build options and dim structures for model inversion
dim.n_theta = 0;
dim.n_phi = 2;
dim.n = 0;
dim.p = 1;
dim.n_t = ntrials;
options.inG = in;
options.dim = dim;
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.DisplayWin = 0;

% Build time series of hidden states and observations
Nmcmc = 16;
nm = length(models);
p = cell(nm,nm,Nmcmc);
o = cell(nm,nm,Nmcmc);
F = zeros(nm,nm,Nmcmc);
for ii=1:Nmcmc
    R = exp(randn(2,ntrials));
    t = exp(randn(2,ntrials));
    u = zeros(4,ntrials);
    u(in.ind.t,:) = t;
    u(in.ind.R,:) = R;
    for i=1:nm
        options.inG.model = models{i};
        [y,x,x0,eta,e] = VBA_simulate (ntrials,[],g_fname,[],phi,u,[],[],options);
        % Call inversion routine
        for j=1:nm
            options.inG.model = models{j};
            [p{i,j,ii},o{i,j,ii}] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);
            F(i,j,ii) = o{i,j,ii}.F;
        end
    end
end


hf = figure('color',[1 1 1]);
ha = axes('parent',hf);
dF = F(:,1,:) - F(:,2,:);
mdF = mean(dF,3);
vdF = var(dF,[],3)./Nmcmc;
plotUncertainTimeSeries(mdF,vdF,[],ha);
set(ha,'xlim',[0,3],'xtick',[1,2],'xticklabels',models)
xlabel(ha,'type of simulated data')
ylabel(ha,['log p(y|',models{1},') - log p(y|',models{2},')'])
box(ha,'off')




