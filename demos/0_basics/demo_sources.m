function [posterior, out, summary] = demo_sources ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_sources ()
% Demo of observation distribution specification
%
% This demo shows how to simulate and invert models having different, and 
% possibly multiple, data distributions
%
% Background
% ~~~~~~~~~~
%
% Let's say you want to measure the ability of your subjects to detect a
% stimulus and you collect three types of observations: the (log-) reaction time
% (continuous, normally distributed), a button press (seen/unseen, binary), 
% and a rating between 1 and 5 (very unlikely <-> very likely,
% multinomial).
% In this simplified example, we assume that all responses are an affine 
% transformation of the stimulus intensity (u). We then map to the observation
% using ad hoc link function (logistic for button press, softmax for
% ratings)
% Although we can fit the respective modalities independently (first part of
% the demo), we can also use the hypothesis that the same generative model
% (ie the same set of parameters) are underlying all observations and
% therefore fit all data at once (second part of the demo) to get a better 
% estimation.
%
% /////////////////////////////////////////////////////////////////////////

%% Define the models
% =========================================================================

% Gaussian observations
% -------------------------------------------------------------------------
function gx = g_source1 (~, phi, ut, ~)
    % affine mapping
    gx = - phi(1) * ut + phi(2);
end

% Binary observations
% -------------------------------------------------------------------------
function gx = g_source2 (~, phi, ut, ~)
    % sigmoid mapping
    gx = VBA_sigmoid(ut, 'slope', phi(1), 'center', phi(2));
end

% Multinomial observations
% -------------------------------------------------------------------------
function gx = g_source3 (~, phi, ut, ~)
    % softmax mapping
    ratings = - 2 : 2;
    gx = exp (- phi(1) * (ut - ratings - phi(2)) .^ 2);
    gx = gx(:) / sum (gx);
end

%% Simulate data 
% =========================================================================

% Experimental design
% -------------------------------------------------------------------------

% generate experimental design
T = 50;
u = randn(1,T);

% Choose basic settings for simulations
phi = [2; 1]; % observation parameters (scale and offset)

% dimension of the models (here the same for all)   
dim.n_phi = 2; % nb of observation parameters


% Simulations
% -------------------------------------------------------------------------

options.sources.type = 0; % flag for gaussian observations (default, can be omited)
sigma = 0.1; % precision of the gaussian noise
y1 = VBA_simulate (T,[],@g_source1,[],phi,u,Inf,sigma,options);

options.sources.type = 1; % flag for Bernoulli observations
y2 = VBA_simulate (T,[],@g_source2,[],phi,u,Inf,[],options);

options.sources.type = 2; % flag for multinomial observations
y3 = VBA_simulate (T,[],@g_source3,[],phi,u,Inf,[],options);

% Inversion
% -------------------------------------------------------------------------
% Here, we will assume that the paramters are NOT shared across observations.
% If we want to compare this to the hypothesis that the same parameters are 
% generating the data, we can simply chain the inversion by carrying over 
% the posterior as a prior for the next inversion. A better way is to
% invert all data at once (see below).

options.sources.type = 0; 
[posterior.s1, out.s1] = VBA_NLStateSpaceModel(y1,u,[],@g_source1,dim,options);

% options.priors = posterior.s1; % <~ uncomment for the shared parameter hypothesis
options.sources.type = 1; 
[posterior.s2, out.s2] = VBA_NLStateSpaceModel(y2,u,[],@g_source2,dim,options);


% options.priors = posterior.s2; % <~ uncomment for the shared parameter hypothesis
options.sources.type = 2; 
[posterior.s3, out.s3] = VBA_NLStateSpaceModel(y3,u,[],@g_source3,dim,options);


%% Combine observations in a joint estimation
% =========================================================================

% concatenate observations
y = [y1; y2; y3];

% define joint observation function: as we need to predicts all lines of
% observations at once, we simply concatenate respective predictions the
% same way we concatenated y1 to y3.
% Note that in general, you don't have to necessarily forward all inputs
% and parameters to each observation subfunctions...
function gx = g_sourceAll (xt, phi, ut, in)
    gx = cat(1, ...
        g_source1 (xt, phi, ut, in), ...
    	g_source2 (xt, phi, ut, in), ...
    	g_source3 (xt, phi, ut, in));
end

options = struct ();
% describe the type of distribution of the different lines of observations
options.sources(1).out = 1; % selecting y1 line in y
options.sources(1).type = 0; % flag for gaussian observations

options.sources(2).out = 2; % selecting y2 line in y
options.sources(2).type = 1; % flag for bernoulli observations

options.sources(3).out = 3 : 7; % selecting y3 lines in y
options.sources(3).type = 2; % flag for multinomial observations

% note that we could generate new multisource data as following:
% y = VBA_simulate (T,[],@g_sourceAll,[],phi,u,Inf,[],options);

% call inversion routine
[posterior.all, out.all] = VBA_NLStateSpaceModel(y,u,[],@g_sourceAll,dim,options);

%% Display results
% =========================================================================

post2str = @ (post) arrayfun (@ (m, v) sprintf('%3.2f (%3.2f)', m, v), post.muPhi, diag(post.SigmaPhi), 'UniformOutput', false);

summary = table ( ...
    post2str (posterior.s1), ...
    post2str (posterior.s2), ...
    post2str (posterior.s3), ...
    mean([posterior.s1.muPhi, posterior.s2.muPhi, posterior.s3.muPhi],2), ...
    post2str (posterior.all), ...
    phi, ...
    'RowNames', {'phi_1', 'phi_2'}, ...
    'VariableNames', {'y1','y2','y3','mean','joint','true'});

fprintf('Posterior parameter estimates:\n')
disp (summary);

p = VBA_softmax ([(out.s1.F + out.s2.F + out.s3.F), out.all.F]);
fprintf ('\nPosterior probability of parameters being shared: %4.3f\n', p(2));


end