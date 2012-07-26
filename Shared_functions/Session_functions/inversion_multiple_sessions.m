function [posterior,out] = inversion_multiple_sessions(in_sessions,...
    y,...
    u,...
    isYout,...
    priors)

% This function performs the inversion of a model for multiple independant
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

dim = in_sessions.dim;

% --- General information for all sessions
in = struct();
in.nsess = in_sessions.dim.n_sess;
in.dim = dim;

% --- Specific information for each session
for i = 1 : in_sessions.dim.n_sess
    
    % Information about the evolution/obsevation function for each session
    in.sess(i).f_fname = in_sessions.f_fname; % the function to be used for each session
    in.sess(i).g_fname = in_sessions.g_fname;
    % Information about indices of parameters, hidden states and output used by
    % each session
    in.sess(i).ind.x = dim.n_ps*(i-1)+1:dim.n_ps*i;
    in.sess(i).ind.gx = dim.p_ps*(i-1)+1:dim.p_ps*i;
    in.sess(i).ind.u = dim.n_data_ps*(i-1)+1:dim.n_data_ps*i;
    
    if isempty(in_sessions.ind.theta)
        in.sess(i).ind.theta = [];
    else in.sess(i).ind.theta = in_sessions.ind.theta(i,:); end
    
    if isempty(in_sessions.ind.phi) 
        in.sess(i).ind.phi = [];
    else in.sess(i).ind.phi = in_sessions.ind.phi(i,:); end
    
    try in.sess(i).inG = in_sessions.inG{i};
    catch; in.sess(i).inG = in_sessions.inG;end
    
    try  in.sess(i).inF = in_sessions.inF{i};
    catch ;in.sess(i).inF = in_sessions.inF;end
    
    
    
end


%%%% Options for inversion
if exist('in_sessions.options','var')
    options = in_sessions.options;
end


options.inF = in;
options.inG = in;

options.priors = priors;

try in_sessions.DisplayWin
    options.DisplayWin = in_sessions.DisplayWin;
catch ; options.DisplayWin = 1; end


options.GnFigs = 0;
options.binomial = in_sessions.binomial; % Dealing with binary data
options.isYout = isYout; % Excluding data points
[posterior,out] = VBA_NLStateSpaceModel(y,u,@f_nsess,@g_nsess,dim,options);

end




