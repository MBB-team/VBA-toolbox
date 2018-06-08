function [h,ep,stat] = bayesian_ttest(x1,x2)
% function [h,ep,stat] = bayesian_ttest(x1,x2)

% bayesian_ttest - performs a statistical test of mean comparison under a
%                  bayesian inferential procedure. The test can be:
%                  1-sample or 2-sample.
%
% Syntax: 
%   [h,ep,stat] = bayesian_ttest(x1) performs a 1-sample test of mean comparison against 0 % for future release
%   [h,ep,stat] = bayesian_ttest(x1,x2) if x2 is a double, performs a 1-sample test of mean comparison against x2 % for future release
%   [h,ep,stat] = bayesian_ttest(x1,x2)    if x2 is a vector, performs a 2-sample test of mean comparison
%
% Inputs:
%    x1 - first sample vector (nsample x 1)
%    x2 - reference mean (1 x 1)
%       - second sample vector (nsample x 1)
%
% Outputs:
%    h - statistical decision to accept the alternative hypothesis (0/1)
%    ep - excedance probability of the alternative hypothesis (from 0 to 1)
%    stat - struct
%       .logBayesFactor - log-bayes factor: log(P(H1|m)/P(H0|m))
%       .posterior - posterior moments of the best hypothesis
%           .mu - mean
%
% Author: Nicolas Borderies
% email: nico.borderies@gmail.com 
% 14-Dec-2017

%% check options
% ----------------------------
if ~exist('x2')
    testType = '1sample';
    muH0 = 0 ;
    error('1-sample test not yet available');
elseif numel(x2)==1
    testType = '1sample';
    muH0 = x2 ;
    error('1-sample test not yet available');
else
    testType  = '2sample';
    x = [x1;x2]';
    mean_x  = VBA_nanmean(x);
end

%% model definitions
% ----------------------------

% inputs/outputs
% ----------------------------
nsamples = [numel(x1) , numel(x2)];
dummy = ones(1,max(nsamples));
samples = nan(2,max(nsamples));
samples(1,1:nsamples(1)) = x1;
samples(2,1:nsamples(2)) = x2;

% formula 
% ----------------------------
g_fname = @g_ttest;

% priors
% ----------------------------
% - 1-sample
%  model : { mu1 = phi_1  }
%  H0: phi_1 = k ; H1: phi_1 ~= k 

% - 2-samples
%  model : { mu1 = phi_1 ; mu2 = phi_1 + phi_2  }
%  H0: phi_2 = 0 ; H1: phi_2 ~= 0 

% observation parameters
priors_1.muPhi = [ VBA_nanmean(x1) ; 0 ];  
priors_1.SigmaPhi = diag([1e3*nanstd(x1),1e3*nanstd(x2)]);      
% observation noise parameters
% (jeffrey's uninformative priors on precision)     
priors_1.a_sigma = [1 1].*1;     
priors_1.b_sigma = [1 1].*1;      
% null-hypothesis
priors_0 = priors_1;
priors_0.SigmaPhi(2,2) = 0;
% store
options.priors = priors_1;       

% dimension
dim.n_phi = 2;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states
dim.p = 2;
dim.n_t = size(samples,2);

% inversion options
%%% source identification
options.sources(1).out=1;
options.sources(1).type=0;
options.sources(2).out=2;
options.sources(2).type=0;
%%% missing data
misssingData = isnan(samples);
options.isYout = zeros(2,max(nsamples));
options.isYout(misssingData) = 1;
samples(misssingData)=0;
%%% others
options.verbose=0;
options.DisplayWin = 0 ;


%% model estimations
% ----------------------------

% H1 full model
[p1,stat1] = VBA_NLStateSpaceModel(samples,dummy,[],g_fname,dim,options);
F1 = stat1.F;

% H0 partial model
p0 = p1;
[F0,p] = VBA_SavageDickey(p1,priors_1,F1,dim,priors_0);
p0.muPhi = p.muPhi;
p0.SigmaPhi = p.SigmaPhi;
p0.a_sigma = p.a_sigma;
p0.b_sigma = p.b_sigma;

%% Statistical inference
% ----------------------------

% test decision 
logBayesFactor = (F1 - F0);
ep = 1./(1+exp(-logBayesFactor));
h = (ep>=0.95);

% posterior estimates
p = VBA_BMA([p0;p1],[F0;F1]);
% posterior.mu = p.muPhi ;
posterior.mu = [ p.muPhi(1) , VBA_sparsifyPrior(p.muPhi(2)) ];

% Formating
stat=struct;
stat.logBayesFactor = logBayesFactor;
stat.posterior = posterior;



end


