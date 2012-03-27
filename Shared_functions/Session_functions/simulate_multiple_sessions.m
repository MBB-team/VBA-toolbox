function  [y,x,x0,eta,e] = simulate_multiple_sessions(in_sessions,...
    u,...
    priors,...
    theta,...
    phi)

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
dim = in_sessions.dim;

for i = 1 : in.nsess
    
        in.sess(i).inG = [];
    in.sess(i).inF = [];
    
    
    in.sess(i).f_fname = in_sessions.f_fname; % the function to be used for each session
    in.sess(i).g_fname = in_sessions.g_fname;

    
    in.sess(i).ind.x = dim.n_ps*(i-1)+1:dim.n_ps*i;
    in.sess(i).ind.gx = dim.p_ps*(i-1)+1:dim.p_ps*i;
    in.sess(i).ind.u = dim.n_data_ps*(i-1)+1:dim.n_data_ps*i;
    
    if isempty(in_sessions.ind.theta)
            in.sess(i).ind.theta = [];
    else
            in.sess(i).ind.theta = in_sessions.ind.theta(i,:);
    end
    
    if isempty(in_sessions.ind.phi)
            in.sess(i).ind.phi = [];
    else
            in.sess(i).ind.phi = in_sessions.ind.phi(i,:);
    end
    
    
end

in.dim_output = dim.p;


%%%% Options for inversion
options.inF = in;
options.inG = in;
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = in_sessions.binomial; % Dealing with binary data
%options.isYout = isYout; % Excluding data points

alpha = Inf;
sigma = Inf;
 [y,x,x0,eta,e] = simulateNLSS(dim.n_t,@f_nsess,@g_nsess,theta,phi,u,alpha,sigma,options,priors.muX0);

end




