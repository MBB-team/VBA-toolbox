function [susceptibility,specificity] = VBA_susceptibility(posterior,out,ind)
% scores the susceptibility of u->y relationships w.r.t. model parameters
% [xiPhi,zetaPhi,xiTheta] = VBA_susceptibility0(out,posterior,ind)
% IN:
%   - out/posterior: the output structures of the VBA inversion
%   - ind: a structure containing the following fields:
%       .phi: indices of observation parameters, whose importance in the
%       u->y relationship one wants to score
%       .theta: idem, for evolution parameters.
% OUT:
%   - xiPhi: nuXnphi relative susceptibility matrix (w.r.t. observation
%   parameters), where nu and nphi are the number of considered inputs and
%   observation parameters, respectively.
%   - xiTheta: nuXntheta relative susceptibility matrix (w.r.t. evolution
%   parameters).

% ========================================================================================
% check parameters
% ========================================================================================

dim = out.dim;

% ________________________________________________________________________________________
% check if models has inputs

if dim.u == 0
    disp('*** Susceptibility analysis: the system does not receive any input!')
    susceptibility = struct('phi',[],'theta',[]);
    specificity    = struct('phi',[],'theta',[]);
    return
end

% ________________________________________________________________________________________
% set default paramters of interest

if ~exist('ind'), ind = struct; end
% check if is DCM
try, isDCM=out.options.inF.fullDCM; catch, isDCM=0; end
% complete option strucutre
if isDCM
    % keep only DCM connections
   idxConnect = 1:(out.options.inF.indself-1);
   ind = check_struct(ind,        ...
        'u'     , 1:dim.u       , ...
        'phi'   , []            , ...
        'theta' , idxConnect      ...
    );  
else
    % all params
    ind = check_struct(ind,       ...
        'u'     , 1:dim.u       , ...
        'phi'   , 1:dim.n_phi   , ...
        'theta' , 1:dim.n_theta   ...
    );  
end

% ========================================================================================
% shortcuts
% ========================================================================================

options = out.options;
f_fname = options.f_fname;
g_fname = options.g_fname;
y = out.y;

% ========================================================================================
% initialization
% ========================================================================================

v0 = score(y,out.suffStat.gx,options) ; 

% ========================================================================================
% Mutilating the models and deriving the output distortion
% ========================================================================================

% ________________________________________________________________________________________
% 1 - switching off inputs
vu = zeros(length(ind.u),1);
for iu=1:length(ind.u)
    mutilated_input_idx = ind.u(iu);
    vu(iu) = mutilation_score(posterior,out,'u',mutilated_input_idx);
end
dvu = v0 - repmat(vu,1,length(ind.phi));

% dealing with parameters
for paramType = {'phi','theta'}
    paramType = paramType{1} ;
    
    if ~isempty(ind.(paramType)) > 0
        
        % compute mutilation score
        vparam  = zeros(length(ind.u),1);
        vuparam = zeros(length(ind.u),length(ind.(paramType)));
        
        % ________________________________________________________________________________
        % 2 - switching off system parameters
        for ip=1:length(ind.(paramType))
            mutilated_param_idx = ind.(paramType)(ip);
            vparam(ip) = mutilation_score(posterior,out,paramType,mutilated_param_idx);
        end
        
        % _______________________________________________________________________________
        % 3 - bilateral mutilations
        for ip=1:length(ind.(paramType))
           for iu=1:length(ind.u)
               mutilated_param_idx = ind.(paramType)(ip);
               mutilated_input_idx = ind.u(iu);
               vuparam(iu,ip) = mutilation_score(posterior,out, ...
                   paramType,mutilated_param_idx,'u',mutilated_input_idx);
           end
        end
        
        % ________________________________________________________________________________
        % compute susceptibility and specificity
        dvuparam = v0 - vuparam;
        % relative susceptibility
        dvphi = v0 - repmat(vparam',length(ind.u),1);
        susceptibility.(paramType) = 1 + dvphi - dvuparam;
        % relative specificity
        specificity.(paramType) = 1 + dvu - dvuparam;
        
    else
        susceptibility.(paramType) = [];
        specificity.(paramType) = [];
    end
    
end

% ========================================================================================
% Additional functions
% ========================================================================================

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% compute the loss of explained variance of a mutilated model
function v = mutilation_score(posterior,out,varargin)

    % default is no mutilation
    mutilation = check_struct(struct,'u',0,'phi',0,'theta',0);
    
    % set mutilation from args
    for iArg = 1:2:numel(varargin)
        mutilation.(varargin{iArg}) = varargin{iArg+1} ;
    end
        
    % compute mutilated model trajectory
    [muy] = laplace_mutilated(posterior,out,mutilation);
    
    % compute score
    v = score(out.y,muy,out.options) ;

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% compute the loss of explained variance of a trajectory
function v = score(y,muy,options)

    % get rid of unused data
    y(options.isYout==1)   = [] ;
    muy(options.isYout==1) = [] ;
    
    % compute loss of explained variance
    dy = vec(y) - vec(muy);
    v = 1 - mean(dy.^2)./var(y) ;
    

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% compute the trajectory moments for a mutilated model
function [muy] = laplace_mutilated(posterior,out,mutilation)
    
    mutilation = check_struct(mutilation,'u',0,'phi',0,'theta',0);

    % define input
    u = out.u;
    if mutilation.u > 0
        u(mutilation.u,:) = 0;
    end
    
    % define parameters
    opt = out.options;
    opt.priors = posterior;
    
    if mutilation.theta > 0
        opt.priors.muTheta(mutilation.theta)      = 0;
        opt.priors.SigmaTheta(mutilation.theta,:) = 0;
        opt.priors.SigmaTheta(:,mutilation.theta) = 0;
    end
    
    if mutilation.phi > 0
        opt.priors.muPhi(mutilation.phi)      = 0;
        opt.priors.SigmaPhi(mutilation.phi,:) = 0;
        opt.priors.SigmaPhi(:,mutilation.phi) = 0;
    end
   
   % compute trajectory moments
   [muy] = VBA_getLaplace(u,out.options.f_fname,out.options.g_fname,out.dim,opt,0);

    
    

