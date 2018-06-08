function [result] = VBA_susceptibility(posterior,out,ind)
% scores the susceptibility of u->y relationships w.r.t. model parameters
% [result] = VBA_susceptibility(posterior,out,ind)
% IN:
%   - out/posterior: the output structures of the VBA inversion
%   - ind: a structure containing the following fields:
%       .u  : indices of inputs of interests
%       .phi: indices of observation parameters, whose importance in the
%       u->y relationship one wants to score
%       .theta: idem, for evolution parameters.
%       .y : indices of observations of interest
% OUT:
%   - susceptibility: nuXnphi and nuxntheta relative susceptibility
%   matrices (w.r.t. observation/evolution parameters), where nu and 
%   nphi/ntehta are the number of considered inputs and parameters, respectively.
%   - susceptibility: nuXntphi qnd nuxntheta relative susceptibility matrices
%   (w.r.t. observation/evolution  parameters).

% ========================================================================================
% check parameters
% ========================================================================================

dim = out.dim;

% ________________________________________________________________________________________
% check if models has inputs

if dim.u == 0
    disp('*** Susceptibility analysis: the system does not receive any input!')
    result.susceptibility = struct('phi',[],'theta',[]);
    result.specificity    = struct('phi',[],'theta',[]);
    return
end

% ________________________________________________________________________________________
% set default parameters of interest

if ~exist('ind'), ind = struct; end

% check if is DCM
try, isDCM=out.options.inF.fullDCM; catch, isDCM=0; end

% complete option strucutre
if isDCM
    % keep only DCM connections
   idxConnect = 1:(out.options.inF.indself-1);
   % idxConnect = setdiff(idxConnect,out.options.inF.indC); % trivial effects (or maybe not)
   if isfield(out.options.inF,'indhself') 
      % behaviour DCM 
      idxResp = [out.options.sources(2:end).out] ;
   else
      % normal DCM
      idxResp = 1:dim.p;
   end
   ind = VBA_check_struct(ind,        ...
        'u'     , 1:dim.u       , ...
        'phi'   , []            , ...
        'theta' , idxConnect    , ...
        'y'     , idxResp         ...
    );  
else
    % all params
    ind = VBA_check_struct(ind,       ...
        'u'     , 1:dim.u       , ...
        'phi'   , 1:dim.n_phi   , ...
        'theta' , 1:dim.n_theta , ...
        'y'     , 1:dim.p         ...
    );  
end

% get rid of fixed parameters
fixedTheta = find(diag(out.options.priors.SigmaTheta)==0);
ind.theta = setdiff(ind.theta,fixedTheta);
fixedPhi = find(diag(out.options.priors.SigmaPhi)==0);
ind.phi = setdiff(ind.phi,fixedPhi);

% get rid of session param
if isfield(out.options, 'multisession')
    ind.u = setdiff(ind.u,dim.u);
end
% ========================================================================================
% shortcuts
% ========================================================================================

%options = out.options;
%f_fname = options.f_fname;
%g_fname = options.g_fname;
%y = out.y;

% ========================================================================================
% initialization
% ========================================================================================



% ========================================================================================
% Mutilating the models and deriving the output distortion
% ========================================================================================

% ________________________________________________________________________________________
% 1 - switching off inputs
vu = zeros(length(ind.u),1);
parfor iu=1:length(ind.u)
    mutilated_input_idx = ind.u(iu);
    vu(iu) = mutilation_score(posterior,out,ind,'u',mutilated_input_idx);
end

% dealing with parameters
for paramType = {'phi','theta'}
    paramType = paramType{1} ;
    
    if ~isempty(ind.(paramType)) > 0
        
        % compute mutilation score
        vparam  = zeros(1,numel(ind.(paramType)));
        % ________________________________________________________________________________
        % 2 - switching off system parameters
        parfor ip=1:length(ind.(paramType))
            mutilated_param_idx = ind.(paramType)(ip);
            vparam(ip) = mutilation_score(posterior,out,ind,paramType,mutilated_param_idx);
        end
        
        % _______________________________________________________________________________
        % 3 - bilateral mutilations 
        vuparam = zeros(numel(ind.u),numel(ind.(paramType)));
        for iu=1:length(ind.u)
           parfor ip=1:length(ind.(paramType))
               mutilated_param_idx = ind.(paramType)(ip);
               mutilated_input_idx = ind.u(iu);
               vuparam_iu(ip) = mutilation_score(posterior,out,ind, ...
               paramType,mutilated_param_idx,'u',mutilated_input_idx);
           end
           vuparam(iu,:) = vuparam_iu ;
        end
        
        % ________________________________________________________________________________
        % compute susceptibility and specificity
        dvuparam = vuparam;
        dvu      = repmat(vu,1,length(ind.(paramType)));
        dvparam  = repmat(vparam,length(ind.u),1);

        % relative susceptibility
        susceptibility.raw.(paramType)  = 1 + (dvparam - dvuparam)      ;
        susceptibility.norm.(paramType) = 1 + (dvparam - dvuparam)./dvu ;
       % relative specificity
        specificity.raw.(paramType)  = 1 + (dvu - dvuparam)      ;
        specificity.norm.(paramType) = 1 + (dvu - dvuparam)./dvparam ;
        
        % interaction scores
        inter =  ((dvparam+dvu)-dvuparam) ; 
        interaction.(paramType) = inter;
        interaction_norm.(paramType)  =  inter./min(dvu,dvparam);
        interaction_normU.(paramType) =  inter./dvu;
        interaction_normP.(paramType) =  inter./dvparam ;
        ii = max(inter,0) ;
        interaction_normC.(paramType) =  ii./ repmat(sum(ii),size(ii,1),1) ;
         
        
    else
        susceptibility.(paramType) = [];
        specificity.(paramType) = [];
        interaction.(paramType) = [];
        interaction_norm.(paramType) = [];
        interaction_normP.(paramType) = [];
        interaction_normU.(paramType) = [];
        interaction_normC.(paramType) = [];
    end
    
end


result.ind = ind ;
result.susceptibility = susceptibility;
result.specificity    = specificity;
result.interaction.raw     = interaction;
result.interaction.norm    = interaction_norm;
result.interaction.normP   = interaction_normP;
result.interaction.normU   = interaction_normU;
result.interaction.normC   = interaction_normC;

end

% ========================================================================================
% Additional functions
% ========================================================================================

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% compute the loss of explained variance of a mutilated model
function v = mutilation_score(posterior,out,ind,varargin)

    % default is no mutilation
    mutilation = VBA_check_struct(struct,'u',0,'phi',0,'theta',0);
    
    % set mutilation from args
    for iArg = 1:2:numel(varargin)
        mutilation.(varargin{iArg}) = varargin{iArg+1} ;
    end
        
    % compute mutilated model trajectory
    [muy] = laplace_mutilated(posterior,out,mutilation);
    
    % compute score
    y = out.y ;
    vfull = score(y,out.suffStat.gx,out.options,ind) ; 
    v = vfull - score(y,muy,out.options,ind) ;
end

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% compute the loss of explained variance of a trajectory
function v = score(y,muy,options,ind)

    % select observation of interest
    y      =   y(ind.y,:);
    muy    = muy(ind.y,:);
    isYout = options.isYout(ind.y,:);
    
    % get rid of unused data
    y(isYout==1)   = [] ; 
    muy(isYout==1) = [] ;
    
    % compute explained variance % ### TODO: multiline data
    dy = VBA_vec(y) - VBA_vec(muy);
    v =  1 - mean(dy.^2)./var(y) ;
end   

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% compute the trajectory moments for a mutilated model
function [muy] = laplace_mutilated(posterior,out,mutilation)
    
    mutilation = VBA_check_struct(mutilation,'u',0,'phi',0,'theta',0);

    mutilation_ratio = 0;
    % define input
    u = out.u;
    if all(mutilation.u > 0)
        u(mutilation.u,:) = mutilation_ratio*u(mutilation.u,:);
    end
    
    % define parameters
    opt = out.options;
    opt.priors = posterior;
    
    if mutilation.theta > 0
        opt.priors.muTheta(mutilation.theta)      = mutilation_ratio*opt.priors.muTheta(mutilation.theta);
        opt.priors.SigmaTheta(mutilation.theta,:) = 0;
        opt.priors.SigmaTheta(:,mutilation.theta) = 0;
    end
    
    if mutilation.phi > 0
        opt.priors.muPhi(mutilation.phi)      = mutilation_ratio*opt.priors.muPhi(mutilation.phi);
        opt.priors.SigmaPhi(mutilation.phi,:) = 0;
        opt.priors.SigmaPhi(:,mutilation.phi) = 0;
    end
   
   % compute trajectory moments
   [muy] = VBA_getLaplace(u,out.options.f_fname,out.options.g_fname,out.dim,opt,0,'skip');
   muy = reshape(muy,size(out.y));
end
    
    

