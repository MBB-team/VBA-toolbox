function [ posterior_s, out_s ] = extract_sessions( sessions, posterior, out )
% Extract inversion information for desired sessions from inversion info for multiple sessions
% INPUT
% - sessions : indices of sessions to extract
% - posterior & out : output of the inversion over multiple sessions from an
% extended model
% OUTPUT
% - posterior_s & out_s : output restricted to the desired sessions.

% NOTE : this corresponds to croping the initial posterior.
% -> This is not a new inversion over choosen sessions
% -> Because we deal with multivariate gaussian variables, this is
% EQUIVALENT TO A MARGINALIZATION of the posterior over parameters not
% involved in the sessions to extract

Nsess = length(sessions); % number of sessions to extract
sess = out.options.inF.sess; % extract information about sessions
dim = out.options.dim; % dimension of inverted model

% -- Extracting indices
ind = out.options.inF.sess(sessions(1)).ind; % first session to create struct
if length(sessions)>1
    for i_sess = sessions(2:end) % remaining sessions
        % merging indices of variables used in the sessions
        ind.x = union(ind.x, sess(i_sess).ind.x);
        ind.theta = union(ind.theta, sess(i_sess).ind.theta);
        ind.phi = union(ind.phi, sess(i_sess).ind.phi);
        ind.gx = union(ind.gx, sess(i_sess).ind.gx);
        ind.u = union(ind.u, sess(i_sess).ind.u);
    end
end

% -- to track changes in indices with model change
% (parameters might be shared by multiple sessions)
phimap = zeros(1,dim.n_phi);
for i = 1 : length(ind.phi);
    phimap(ind.phi(i))=i;
end
thetamap = zeros(1,dim.n_theta);
for i = 1 : length(ind.theta);
    thetamap(ind.theta(i))=i;
end

% -- Calculating new dimensions
dim.n = length(ind.x);
dim.n_theta = length(ind.theta);
dim.n_phi = length(ind.phi);
dim.p = length(ind.gx);
dim.u = length(ind.u);
dim.n_sess = Nsess;

out_s = out;
out_s.options.dim = dim;
out_s.dim = dim;

% -- Calculating new posterior
posterior_s = posterior;
posterior_s.muPhi = posterior.muPhi(ind.phi);
posterior_s.muTheta = posterior.muTheta(ind.theta);
posterior_s.SigmaPhi = posterior.SigmaPhi(ind.phi,ind.phi);
posterior_s.SigmaTheta = posterior.SigmaTheta(ind.theta,ind.theta);

posterior_s.muX0 = posterior.muX0(ind.x);
posterior_s.SigmaX0 = posterior.SigmaX0(ind.x,ind.x);
posterior_s.muX = posterior.muX(ind.x,:);
for i = 1 : dim.n_t
    posterior_s.SigmaX.current{i} = posterior.SigmaX.current{i}(ind.x,ind.x);
end

% -- Calculating new priors
priors_s = out.options.priors;
priors_s.muPhi = out.options.priors.muPhi(ind.phi);
priors_s.muTheta = out.options.priors.muTheta(ind.theta);
priors_s.SigmaPhi = out.options.priors.SigmaPhi(ind.phi,ind.phi);
priors_s.SigmaTheta = out.options.priors.SigmaTheta(ind.theta,ind.theta);

priors_s.muX0 = out.options.priors.muX0(ind.x);
priors_s.SigmaX0 = out.options.priors.SigmaX0(ind.x,ind.x);
out_s.options.priors = priors_s;

% -- Calculating new suffStat
out_s.suffStat = out.suffStat;
out_s.suffStat.gx = out.suffStat.gx(ind.gx,:);
out_s.suffStat.vy = out.suffStat.vy(ind.gx,:);
out_s.suffStat.dphi = out.suffStat.dphi(ind.phi);
out_s.suffStat.dtheta = out.suffStat.dtheta(ind.theta);
out_s.suffStat.dy = out.suffStat.dy(ind.gx,:);
for i = 1 : dim.n_t
    out_s.suffStat.SigmaX{i} = out.suffStat.SigmaX{i}(ind.x,ind.x);
end

% -- Declaring new options

out_s.options.isYout = out.options.isYout(ind.gx,:);

% -- Calculating new data
out_s.y = out.y(ind.gx,:);

% -- Declaring new functions
if Nsess == 1
    out_s.options.f_fname = sess(sessions(1)).f_fname;
    out_s.options.g_fname = sess(sessions(1)).g_fname;
else % Nsess > 1 => create extended model for multiple sessions
    
    dim = cell(1,Nsess); % dimension of models for each individual session 
    in_sessions = struct();
    in_sessions.n_sess = Nsess; % number of sessions
    in_sessions.dim_e = out_s.dim; % specify extended model's parameter space
    
    for i =  1 : Nsess
        in_sessions.f_fname{i} =  sess(sessions(i)).f_fname;
        in_sessions.g_fname{i} =  sess(sessions(i)).g_fname;
        dim = sess(sessions(i)).dim; % dimension of models for each individual session      
        in_sessions.inF{i} = sess(sessions(i)).inF;
        in_sessions.inG{i} = sess(sessions(i)).inG;
        in_sessions.ind.theta(i,:) =  thetamap(sess(sessions(i)).ind.theta); % specify parameter use for each session
        in_sessions.ind.phi(i,:) =    phimap(sess(sessions(i)).ind.phi); % specify parameter use for each session
    end
    [ f_fname_e,g_fname_e,dim_e,options_e ] = makeExtendedModel(dim,out_s.options,in_sessions); % building new model for multiple sessions
    out_s.options = options_e;
    
end

