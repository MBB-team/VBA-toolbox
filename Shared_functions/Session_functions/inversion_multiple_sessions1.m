function [posterior,out] = inversion_multiple_sessions1(in_sessions,...
    y,...
    u,...
    isYout,...
    priors)

% - in_sessions : information about session
%       - .f_fname : evolution function
%       - .g_fname : observation function
%       - .dim : dimensions of the variables of the model
%       - .ind.theta : indices of the variable theta used for each session
%       - .ind.phi : indices of the variable phi used for each session
% - y : ((dim_output*Nsession)*Ntrials) Observed behavior for inversion
% - u : ((dim_data*Nsession)*Ntrials) Experimenter data needed for model inversion
% - isYout : ((dim_output*Nsession)*Ntrials) behavioral data not to be
%           considered for inversion (1=out,0=in)
% - priors: a structure containing the parameters of the prior pdf, ie:
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


in = struct();
in.nsess = in_sessions.dim.n_sess;
in.f_fname = in_sessions.f_fname;
in.g_fname = in_sessions.g_fname;
dim = in_sessions.dim;
in.dim = dim;
in.ind_theta = in_sessions.ind.theta;
in.ind_phi = in_sessions.ind.phi;
in.dim_output = dim.p;
in.dim_output_ps = dim.p_ps;


%%%% Options for inversion
options.inF = in;
options.inG = in;
options.priors = priors;
try in_sessions.DisplayWin
     options.DisplayWin = in_sessions.DisplayWin;
catch
     options.DisplayWin = 1;
end
options.GnFigs = 0;
options.binomial = in_sessions.binomial; % Dealing with binary data
options.isYout = isYout; % Excluding data points

[posterior,out] = VBA_NLStateSpaceModel(y,u,@f_nsess1,@g_nsess1,dim,options);

end




