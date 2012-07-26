close all 
clear all

% Wolpert Körding 2004 : Bayesian Integration in Sensorimotor learning
% This scripts simulates optimal integration of prior knowledge and noisy
% data as in the aforementioned article
% - A simulation is performed with some prior knowledge about the shift.
% - Inversion from generated data is performed : parameters of the learned gaussian
% prior on shift is infered from behavior

% Remark : there is no learning here. Prior was learned from prior
% training with feedback.

% Simulating data
Ntrials = 200; % Number of time trials
shifts = 1+0.5*randn(1,Ntrials); %sampling shifts
i_cues = randi(4,1,Ntrials); %index of cues
std_cues = zeros(1,Ntrials); %different std of sensed cue
std_cues(find(i_cues == 1)) = 1;
std_cues(find(i_cues == 2)) = 2;
std_cues(find(i_cues == 3)) = 1e5;
std_cues(find(i_cues == 4)) = 0;

% similating model behavior
in.shifts = shifts;
in.std_cues = std_cues;
phi = [1;log(0.5)];
[gx] = g_wk2004([],phi,[],in) + randn(1,Ntrials)*0.05; % estimated shift + gaussian noise

% plotting results
figure
hold on
col = 'rgbc';
xmin = min(shifts);ymin = min(gx);
xmax = max(shifts);ymax = max(gx);
xt = xmin :0.1:xmax;
for i = 1 : 4
    I = find(i_cues == i);
    p = polyfit(shifts(I),gx(I),1);
    plot(shifts(I),gx(I),['x',col(i)])
end
title('Real shifts against estimated shifts for each cue condition')
legend('cue 1','cue 2','cue 3', 'cue 4')
for i = 1 : 4
    I = find(i_cues == i);
    p = polyfit(shifts(I),gx(I),1);
    plot([xmin,xmax],p(2)+p(1)*[xmin,xmax],col(i))
end
axis([xmin,xmax,ymin,ymax])

%%
% Inverting the model
g_fname = @g_wk2004; 
dim = struct('n',0,'n_theta',0,'n_phi',2,'p',Ntrials,'n_t',1);
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.SigmaPhi =1e2*eye(dim.n_phi);
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;
% Options for inversion
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 0; % Dealing with continuous data
options.dim = dim;
options.verbose = 0;
options.inG = in;
y = gx';

[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);  % Inversion function
% Estimated priors 
disp('-------- Infered prior about shift')
disp(['Estimated mean of prior on shift : ',num2str(posterior.muPhi(1))])
disp(['Estimated std of prior on shift : ',num2str(exp(posterior.muPhi(2)))])

displayResults(posterior,out,y,[],[],[],phi,Inf,0)

