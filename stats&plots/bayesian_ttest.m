function [h,ep,stat] = bayesian_ttest(x1,x2)
% bayesian_ttest - performs a statistical test of mean comparison under a
%                  bayesian inferential procedure. The test can be:
%                  1-sample or 2-sample.
%
% Syntax: 
%   [h,ep,stat] = bayesian_ttest(x1) performs a 1-sample test of mean comparison against   0 % for future release
%   [h,ep,stat] = bayesian_ttest(x1,x2) if x2 is a double, performs a 1-sample test of mean comparison against x2 % for future release
%   [h,ep,stat] = bayesian_ttest(x1,x2)    if x2 is a vector, performs a 2-sample test of mean comparison
%
% Inputs:
%    x1 - first sample vector (nsample x 1)
%    x2 - reference mean (1 x 1)
%       - second sample vector (nsample x 1)
%       
%
% Outputs:
%    h - statistical decision to accept the alternative hypothesis (0/1)
%    ep - excedance probability of the alternative hypothesis (from 0 to 1)
%    stat - struct
%       .lppRatio - log-posterior probability ratio (H1/H0)
%       .posterior - posterior moments of the best hypothesis
%           .mu - mean
%           .sigma - standard deviation
%
%
% Example: 
%
%
% Requirements:
%   Script: 
%   Subfunctions: sem.m
%   Data-files: 
%   Matlab-version:
%   Matlab-toolbox: VBA-toolbox
%
% See also: 

% Author: Nicolas Borderies
% email address: nico.borderies@gmail.com 
% February 2017; Last revision: 

%% -------------------------------------

% check options
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
    muH0 = VBA_nanmean(x);
    muH1 = [VBA_nanmean(x1),VBA_nanmean(x2)];
    sigmaH0 = sem(x);
    sigmaH1 = [sem(x1),sem(x2)];

end

% inputs/outputs
nsamples = [numel(x1) , numel(x2)];
dummy = ones(1,max(nsamples));
samples = nan(2,max(nsamples));
samples(1,1:nsamples(1)) = x1;
samples(2,1:nsamples(2)) = x2;

% formula 
g_fname = @g_ttest;

% priors
priors.muPhi = repmat(muH0,2,1);        % H0: mu1=mu2 
% priors.SigmaPhi = sigmaH0*ones(2);      % H0: mu1=mu2 
priors.SigmaPhi = sigmaH0*eye(2);       % H0: mu1=mu2 
priors.a_sigma = [1 1];             % Jeffrey's prior
priors.b_sigma = [1 1];                % Jeffrey's prior
options.priors = priors;       

% dimension
dim.n_phi = 2;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states

% inversion options
options.sources(1).out=1;
options.sources(1).type=0;
options.sources(2).out=2;
options.sources(2).type=0;
options.verbose=0;
misssingData = isnan(samples);
options.isYout = zeros(2,max(nsamples));
options.isYout(misssingData) = 1;
samples(misssingData)=0;
options.DisplayWin = 0 ;


% Estimate posterior densities of both hypothesis
%%% H0:null hypothesis
inG.hypothesis = 'null'; % H0: mu1=mu2 
options.inG = inG;
options.priors = priors;
[p0,stat0] = VBA_NLStateSpaceModel(samples,dummy,[],g_fname,dim,options);

%%%% H1:null hypothesis
inG.hypothesis = 'alternative'; % H1: mu1~=mu2 
options.inG = inG;
options.priors = priors;
[p1,stat1] = VBA_NLStateSpaceModel(samples,dummy,[],g_fname,dim,options);

% Perform statistical comparison
lppRatio = (stat1.F - stat0.F);
ep = 1./(1+exp(-lppRatio));
h = (ep>=0.95);
if lppRatio>=0
    posterior.mu = p1.muPhi;
    posterior.sigma = [p1.SigmaPhi(1,1);p1.SigmaPhi(2,2)];
else
    posterior.mu = p0.muPhi(1);
    posterior.sigma = p0.SigmaPhi(1,1);
end

% Formating
stat=struct;
stat.lppRatio = lppRatio;
stat.posterior = posterior;



end

function [ stdErrMean ] = sem(normal_sample,dim)
%sem: standard error of the mean (for sample following a normal law)
%     input : - sample observed saved into a vector
%             - dimension of the vector indicating index of observation
%
% default argument
if nargin == 1
    dim = find(size(normal_sample) == max(size(normal_sample)),1,'first' );
end
stdErrMean = VBA_nanstd(normal_sample,0,dim)/sqrt(size(normal_sample,dim));


end
