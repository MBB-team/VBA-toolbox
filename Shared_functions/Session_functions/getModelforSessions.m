function [ f_fname_e,g_fname_e,dim_e,options_e,fb_e ] = getModelforSessions( f_fname,g_fname,dim,options,in_sessions)

% This function generates a model to invert multiple sessions.

%{
 This function performs the inversion of a model for multiple independant
% sessions.
% All sessions are independent results of the same model (called base model).
% Each session has different hidden states but can share parameters.
%
% A bigger model encapsulating multiple base models is constructed (called
% extended model). This model generates all sessions independently.
% The extended model calls functions of the base model for subsets of its
% hidden states and parameters.
% The subsets of parameters called for each session are declared as input.
% The priors on parameters and hidden states are declared as input

% - in_sessions : information about session
%       - .f_fname : evolution function for a single session (base model)
%       - .g_fname : observation function for a single function (base
%       model)
%       - .dim : dimensions of the variables of the model for a single
%       session (base model)
%       - .ind.theta : indices of the variable theta used for each session
%       - .ind.phi : indices of the variable phi used for each session
% - y : ((dim_output*Nsession)*Ntrials) Output of all sessions concatenated
% - u : ((dim_data*Nsession)*Ntrials) Experimenter data of all sessions
%       concatenated
% - isYout : ((dim_output*Nsession)*Ntrials) behavioral data not to be
%           considered for inversion (1=out,0=in), concatenated for all
%           sessions
% - priors: a structure containing the parameters of the prior pdf of the
% extended model :
%       .muPhi: a n_phix1 vector containing the prior mean of Phi, the
%       observation parameters
%       .muTheta: a n_thetax1 vector containing the prior mean of Theta,
%       the evolution parameters
%       .muX0: a nx1 vector containing the prior mean of the hidden-states
%       initial condition
%       .SigmaPhi: n_phixn_phi prior covariance matrix of Phi
%       .SigmaTheta: n_thetaxn_theta prior covariance matrix of Theta
%       .SigmaX0: nxn prior covariance matrix of X0
%       .a_sigma / .b_sigma: the shape and scale parameters of the prior
%       Gamma pdf upon the measurement noise precision
%       .a_alpha / .b_alpha: the shape and scale parameters of the prior
%       Gamma pdf upon the stochastic innovations precision

%}


n_sess = in_sessions.n_sess;

%---------------------------------------------------
%-- Dimensions of the extended model


dim_e = struct();
dim_e.p = dim.p*n_sess; % output
dim_e.n_t = dim.n_t; % number of trials (unchanged)
dim_e.u = dim.u*n_sess; % number of trials (unchanged)
dim_e.n = dim.n*n_sess; % hidden states
dim_e.n_sess = n_sess;

%-- Parameters
%- theta
try
    dim_e.n_theta = in_sessions.dim_e.n_theta;
catch
    dim_e.n_theta = dim.n_theta*n_sess;
    in_sessions.ind.theta = reshape(1:dim_e.n_theta,n_sess,dim.n_theta);
end
%- phi
try
    dim_e.n_phi = in_sessions.dim_e.n_phi;
catch
    dim_e.n_phi = dim.n_phi*n_sess;
    in_sessions.ind.phi = reshape(1:dim_e.n_phi,n_sess,dim.n_phi);
end

%------------------------------------------------------
%-- Filling session structure

% --- General information for all sessions
in = struct();
in.nsess = in_sessions.n_sess;


% --- Specific information for each session
for i = 1 : n_sess
    
    % Information about the evolution/obsevation function for each session
    try in.sess(i).f_fname = in_sessions.f_fname{i};
    catch; in.sess(i).f_fname = in_sessions.f_fname;end
    try  in.sess(i).g_fname = in_sessions.g_fname{i};
    catch ;in.sess(i).g_fname = in_sessions.g_fname;end
    
    
    % Information about the evolution/obsevation function for each session
    in.sess(i).f_fname = in_sessions.f_fname; % the function to be used for each session
    in.sess(i).g_fname = in_sessions.g_fname;
    
    % Information about indices of parameters, hidden states and output used by
    % each session
    in.sess(i).ind.x = dim.n*(i-1)+1:dim.n*i;
    in.sess(i).ind.gx = dim.p*(i-1)+1:dim.p*i;
    in.sess(i).ind.u = dim.u*(i-1)+1:dim.u*i;
    
    % Information about the evolution/obsevation paramaters for each session
    if isempty(in_sessions.ind.theta)
        in.sess(i).ind.theta = [];
    else in.sess(i).ind.theta = in_sessions.ind.theta(i,:); end
    if isempty(in_sessions.ind.phi)
        in.sess(i).ind.phi = [];
    else in.sess(i).ind.phi = in_sessions.ind.phi(i,:); end
    
    % Information about the evolution/obsevation extra input for each session
    try in.sess(i).inG = in_sessions.inG{i};
    catch; in.sess(i).inG = in_sessions.inG;end
    try  in.sess(i).inF = in_sessions.inF{i};
    catch ;in.sess(i).inF = in_sessions.inF;end
    
end

%------------------------------------------------------
%---- Options

options_e = struct();
options_e.inF = in;
options_e.inG = in;
options_e.inG.dim = dim_e;

options_e.GnFigs = 0;
options_e.binomial = in_sessions.binomial; % Dealing with binary data
try  options_e.DisplayWin = in_sessions.DisplayWin;
catch ; options_e.DisplayWin = 1; end
try  options_e.binomial = in_sessions.binomial;
catch ; disp('You haven''t specified the type of the data in the options : options.binomial')
end

%---- Function handles

f_fname_e = @f_nsess;
g_fname_e = @g_nsess;


%---- Handle for simulation
fb_e = [];
try
    fb = options.fb;
    fb_e = struct('h_fname',@h_nsess,...
        'nsess',n_sess);
    for i = 1 : n_sess
        fb_e.inH.sess(i).h_fname = fb.h_fname;
        fb_e.inH.sess(i).indy = fb.indy + dim.p*(i-1);
        fb_e.inH.sess(i).indgx = dim.p*(i-1)+1:dim.p*i;

       % fb_e.inH.sess(i).indy = fb.indy + dim.p*(i-1);
        fb_e.inH.sess(i).indfb =  fb.indfb+  dim.u*(i-1);
        fb_e.inH.sess(i).inH =  fb.inH;

    end
    fb_e.inH.dim_fb = length(fb.indfb)*n_sess;
    fb_e.inH.dim_singlefb = length(fb.indfb);

    fb_e.inH.nsess = n_sess;
    fb_e.indy = fb.indy + dim.u*([1 : n_sess]-1); % this is where to store y in u
    fb_e.indfb =  repmat(fb.indfb,n_sess,1) +  dim.u*([0:n_sess-1]')*(ones(size(fb.indfb)))%-1%([1 : n_sess]-1);

    
    options_e.fb = fb_e;
catch
end
 

end






