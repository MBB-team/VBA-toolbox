function [posterior_sub,out_sub,posterior_group,out_group] = VBA_MFX(y,u,f_fname,g_fname,dim,options)
% VB treatment of mixed-effects analysis
% function [posterior,out] = VBA_MFX(y,u,f_fname,g_fname,dim,options)
% This function approaches model inversion from an empirical Bayes
% perspective, whereby within-subject priors are iteratively refined and
% matched to the inferred parent population distribution.
%  Note: all subjects must use the same model
% IN:
%   - y: nsx1 cell array of observations, where ns is the number of
%   subjects in the group
%   - u:  nsx1 cell array of inputs
%   - f_fname/g_fname: evolution/observation function handles
%   - dim: structure containing the model dimensions.
%   - options: nsx1 cell array of options structure. Note: if specified
%   here, the priors on observation and evolution parameters (as well as
%   initial conditions) are useless, since they are replaced by empirical
%   Bayes priors. Priors on precision hyperparameters however, are not
%   treated as random effects drawn from a parent population distribution.
%   In turn, MFX analysis does not update their moments using group-level
%   information...
%   - priors_group: structure containing the prior sufficient statistics on
%   the moments of the parent population distributions (for observation and
%   evolution parameters, as well as for initial conditions, if
%   applicable). See p_group subfields below.
%   - options_group: Structure containing options for MFX. Fields:
%     .TolFun - Minimum change in the free energy, default is 2e-2
%     .MaxIter - Maximum number of iterations, default is 16
%     .DisplayWin - Do we want a graphical output?
%     .verbose - Do we want verbose text output?
%
% OUT:
%   - p_sub/o_sub: nsx1 cell arrays containng the VBA outputs of the
%   within-subject model inversions.
%   - p_group: structure containing the sufficient statistics of the
%   posterior over the moments of the parent population distribution. Its
%   subfields are:
%       .muPhi/SigmaPhi: VB sufficient statistics (first 2 moments) of the
%       Gaussian posterior pdf over the population mean of observation
%       parameters.
%       .muTheta/SigmaTheta: [id] for evolution parameters.
%       .muX0/SigmaX0: [id] for initial conditions.
%       .a_vPhi/b_vPhi: VB sufficient statistics (scale and shape
%       parameters) of the Gamma posterior pdf over the population
%       precision of observation parameters. NB: a_vPhi and b_vPhi have the
%       same dimension than muPhi!
%       .a_vTheta/b_vTheta: [id] for evolution parameters.
%       .a_vX0/b_vX0: [id] for initial conditions.
%   - o_group: output structure of the VBA_MFX approach. In particular, it
%   contains the following subfields:
%       .F: a vector of free energies (across VB iterations). Its last
%       entry (F(end)) provides the free energy lower bound to the MFX
%       model.
%       .it: the final number of VB iterations
%       .date: date vector for track keeping
%       .initVBA: a structure containing the VBA outputs of the
%       within-subject model inversions, without MFX-type priors
%       (initialization).


%% Check parameters
% =========================================================================
if ~ exist('options','var')
    options = struct ();
end
options = VBA_check_struct (options, VBA_defaultOptions ());

%% Shortcuts
% =========================================================================
% Number of subjects
nS = length(y); 

% expand inputs to an array if necessary
if isnumeric (u)
    temp = cell (1, nS);
    temp(:) = {u};
    u = temp;
end

% expand number of observations if necessary
if numel (dim.n_t) == 1
    n_t = repmat (dim.n_t, 1, nS);
else
    n_t = dim.n_t;
end

%% Check priors
% =========================================================================
% Default priors are used if priors are not explicitly provided through the
% priors_group structure. This means Gaussian(0,1) priors for the
% population mean of observation/evolution parameters and initial
% conditions, and Gamma(1,1) for the corresponding population precisions.
if ~ isfield(options,'priors')
    options.priors = struct ();
end
priors_group = VBA_check_struct (options.priors, VBA_defaultMFXPriors (dim));

%% Initialization
% =========================================================================
% Here, we simply initialize the posterior on the population's mean and
% precision over observation/evolution parameters and initial conditions
% using their prior.

% display
VBA_disp ('VBA treatment of MFX analysis: initialization...', options)

% indices of fixed, random, and locked effects
ind.phi_ffx = find(isFixedEffect(priors_group.a_vPhi,priors_group.b_vPhi));
ind.phi_rfx = find(~ isFixedEffect(priors_group.a_vPhi,priors_group.b_vPhi));
ind.phi_in = find(diag(priors_group.SigmaPhi)~=0);

ind.theta_ffx = find(isFixedEffect(priors_group.a_vTheta,priors_group.b_vTheta));
ind.theta_rfx = find(~ isFixedEffect(priors_group.a_vPhi,priors_group.b_vPhi));
ind.theta_in = find(diag(priors_group.SigmaTheta)~=0);

ind.x0_ffx = find(isFixedEffect(priors_group.a_vX0,priors_group.b_vX0));
ind.x0_rfx = find(~ isFixedEffect(priors_group.a_vX0,priors_group.b_vX0));
ind.x0_in = find(diag(priors_group.SigmaX0)~=0);

out_group.ind = ind;

% options
options.dim = dim;
options.dim.ns = nS;

%% Evaluate within-subject free energies under the prior
% =========================================================================

% loop over subjects
for i = 1 : nS
    
    options_subject{i}.priors = getSubjectPriors(priors_group, ind, nS);

    % VBA model inversion
    options_subject{i}.DisplayWin = 0;
    options_subject{i}.verbose = 0;
    options_subject{i}.MaxIter = 0;
    
    % subject's number of observation
    dim.n_t = n_t(i);  
    
    % dummy inversion
    [posterior_sub(i),out_sub(i)] = VBA_NLStateSpaceModel(y{i},u{i},f_fname,g_fname,dim,options_subject{i});
    
    % store options for future inversions
    options_subject{i} = out_sub(i).options;
    options_subject{i}.MaxIter = 32;

end

posterior_group = priors_group;
out_group.options = options;

F(1) = MFX_F (posterior_sub, out_sub, posterior_group, priors_group, dim, ind);

out_group.F = F;
out_group.it = 0;
out_group.options = options;

[out_group.options] = VBA_displayMFX (posterior_sub, out_sub, posterior_group, out_group, 0, 'off');

%% Iterate VB until convergence
% =========================================================================
% We now update the within-subject effects as well as respective population
% moments according to the mean-field VB scheme. This effectively
% iteratively replaces the priors over within-subject effects by the VB
% estimate of the group mean and precision. The free energy of the ensuing
% MFX procedure is computed for tracking algorithmic convergence.

fprintf (1, ['Main VB inversion...']);

for it = 1 : options.MaxIter
    
    % perform within-subject model inversions
    for i = 1 : nS
        
        % display
        try
            set(out_group.options.display.ho,'string',['VB iteration #',num2str(it),': within-subject model inversions (',num2str(floor(100*(i-1)/nS)),'%)'])
        end
        
        % re-define within-subject priors
        options_subject{i}.priors = updateSubjectPriors (options_subject{i}.priors, posterior_group, ind);

        % bypass VBA initialization
        in{i}.posterior = posterior_sub(i);
        in{i}.out = out_sub(i);
        in{i}.out.options = options_subject{i};
        
        % VBA model inversion
        [posterior_sub(i), out_sub(i)] = VBA_NLStateSpaceModel(y{i},u{i},f_fname,g_fname,dim,options_subject{i},in{i});
          
    end
    
    % display
    try
        set(out_group.options.display.ho,'string',['MFX: updating moments of parent distribution...'])
    end
    
    % update moments of the parent population distribution
    if dim.n_phi > 0
        [posterior_group.muPhi,posterior_group.SigmaPhi,posterior_group.a_vPhi,posterior_group.b_vPhi] = ...
            MFX_VBupdate ( ...
                priors_group.muPhi, ...
                priors_group.SigmaPhi, ...
                {posterior_sub.muPhi}, ...
                {posterior_sub.SigmaPhi}, ...
                posterior_group.a_vPhi, ...
                posterior_group.b_vPhi, ...
                priors_group.a_vPhi, ...
                priors_group.b_vPhi, ...
                ind.phi_ffx, ...
                ind.phi_in ...
            );
    end
    if dim.n_theta > 0
        [posterior_group.muTheta,posterior_group.SigmaTheta,posterior_group.a_vTheta,posterior_group.b_vTheta] = ...
            MFX_VBupdate ( ...
                priors_group.muTheta, ...
                priors_group.SigmaTheta, ...
                {posterior_sub.muTheta}, ...
                {posterior_sub.SigmaTheta}, ...
                posterior_group.a_vTheta, ...
                posterior_group.b_vTheta, ...
                priors_group.a_vTheta, ...
                priors_group.b_vTheta, ...
                ind.theta_ffx, ...
                ind.theta_in ...
            );
    end
    if dim.n > 0
        [posterior_group.muX0,posterior_group.SigmaX0,posterior_group.a_vX0,posterior_group.b_vX0] = ...
            MFX_VBupdate ( ...
                priors_group.muX0, ...
                priors_group.SigmaX0, ...
                {posterior_sub.muX0}, ...
                {posterior_sub.SigmaX0}, ...
                posterior_group.a_vX0, ...
                posterior_group.b_vX0, ...
                priors_group.a_vX0, ...
                priors_group.b_vX0, ...
                ind.x0_ffx, ...
                ind.x0_in ...
            );
    end
    
    % free energy
    F(it + 1) = MFX_F(posterior_sub,out_sub,posterior_group,priors_group,dim,ind);
    out_group.F = F;

    % store initial within-subject VBA model inversion
    if it == 1
        out_group.initVBA.p_sub = posterior_sub;
        out_group.initVBA.o_sub = out_sub;
    end
    
    % display
    [out_group.options] = VBA_displayMFX(posterior_sub,out_sub,posterior_group,out_group);
    
    % if convergence, no need to continue iterating
    dF = F(it+1) - F(it);
    if abs(dF) <= options.TolFun 
        break;
    end
    
end

% Wrapping up
% =========================================================================

    % keep track of iterations
    out_group.it = it;

    
fprintf([' done.','\n'])
out_group.date = clock;
out_group.options.sources = out_sub(1).options.sources;
[out_sub.diagnostics] = deal(NaN);
for i = 1 : nS
    out_group.within_fit.F(i) = out_sub(i).F(end);
    out_group.within_fit.R2(i,:) = out_sub(i).fit.R2;
    out_group.within_fit.LLH0(i) = VBA_LMEH0(out_sub(i).y,out_sub(i).options);
    [tmp,out_sub(i)] = VBA_getDiagnostics(posterior_sub(i),out_sub(i));
end
[out_group.options] = VBA_displayMFX(posterior_sub,out_sub,posterior_group,out_group);
try
    if floor(out_group.dt./60) == 0
        timeString = [num2str(floor(out_group.dt)),' sec'];
    else
        timeString = [num2str(floor(out_group.dt./60)),' min'];
    end
    set(out_group.options.display.ho,'string',['VB treatment of MFX analysis complete (took ~',timeString,').'])
end
try
    str = VBA_summaryMFX(out_group);
    VBA_disp(str,opt)
end
out_group.options.display = [];

end

%% ########################################################################
%  Subroutines
% #########################################################################

%% Variational step of the group level statistics
% =========================================================================
function [m, V, a, b] = MFX_VBupdate (mu_0, V_0, mu_s, Vs, a, b, a0, b0, iFfx, iIn)

    % initialisation
    % ---------------------------------------------------------------------
    % number of subjects
    n_s = numel (mu_s);
    % number of parameters
    n = size (mu_0, 1);
    % inverse prior precision
    iV_0 = VBA_inv (V_0);

    % Locked effects (no update)
    % ---------------------------------------------------------------------
    m = mu_0;
    V = iV_0;
    
    % Random effects
    % ---------------------------------------------------------------------
    % parameter indices
    iRfx = setdiff (1 : n, iFfx);
    iRfx = intersect (iRfx, iIn);

    if ~ isempty (iRfx)
        
        % posterior over group mean
        iQ = diag (a(iRfx) ./ b(iRfx));
        sm = sum ([mu_s{:}], 2); 

        V(iRfx, iRfx) = VBA_inv (iV_0(iRfx, iRfx) + n_s * iQ);
        m(iRfx) = V(iRfx, iRfx) * (iV_0(iRfx, iRfx) * mu_0(iRfx) + iQ * sm(iRfx));

        % posterior of group variance
        sv = 0;
        for i = 1 : n_s
            sv = sv ...
                + (mu_s{i}(iRfx) - mu_0(iRfx)) .^ 2 ...
                + diag (Vs{i}(iRfx,iRfx));
        end

        a(iRfx) = a0(iRfx) + 0.5 * n_s;
        b(iRfx) = b0(iRfx) + 0.5 * (sv + n_s * diag (V(iRfx, iRfx))); % do not index sv
    end
    
    % Fixed effects
    % ---------------------------------------------------------------------
    % parameter indices
    iFfx = intersect (iFfx, iIn);

    if ~ isempty (iFfx)

        % posterior over group mean
        wsm = 0;
        sum_iVs = 0;

        for i = 1 : n_s
            iVs = VBA_inv (Vs{i});
            wsm = wsm + iVs * mu_s{i};
            sum_iVs = sum_iVs + iVs;
        end
        tmp = VBA_inv (sum_iVs);
        V(iFfx, iFfx) = tmp(iFfx, iFfx);
        m(iFfx) = V(iFfx, iFfx) * wsm(iFfx);
        
        % posterior of group variance (not updated, by definition)
        a(iFfx) = a0(iFfx);
        b(iFfx) = b0(iFfx);
    end
end

%% Free energy of the full hierarchical model
% =========================================================================
function F = MFX_F (posterior_sub, out_sub, posterior_group, priors_group, dim, ind)

    ns = length (posterior_sub);

    % within subject free energy
    F = sum ([out_sub.F]);

    % group level correction terms
    if dim.n_phi > 0
        F = F + FreeEnergy_var ( ...
                ns, ...
                posterior_group.muPhi, posterior_group.SigmaPhi,...
                priors_group.muPhi, priors_group.SigmaPhi,...
                posterior_group.a_vPhi,posterior_group.b_vPhi,...
                priors_group.a_vPhi, priors_group.b_vPhi,...
                ind.phi_ffx, ind.phi_in ...
            );
    end
    if dim.n_theta > 0
        F = F + FreeEnergy_var ( ...
                    ns, ...
                    posterior_group.muTheta, posterior_group.SigmaTheta,...
                    priors_group.muTheta, priors_group.SigmaTheta,...
                    posterior_group.a_vTheta, posterior_group.b_vTheta,...
                    priors_group.a_vTheta, priors_group.b_vTheta,...
                    ind.theta_ffx, ind.theta_in ...
                );
    end
    if dim.n > 0
        F = F + FreeEnergy_var ( ...
                    ns,...
                    posterior_group.muX0, posterior_group.SigmaX0,...
                    priors_group.muX0, priors_group.SigmaX0,...
                    posterior_group.a_vX0, posterior_group.b_vX0,...
                    priors_group.a_vX0, priors_group.b_vX0,...
                    ind.x0_ffx, ind.x0_in ...
                );
    end
end

% variable-specific free energy group-level correction term
% =========================================================================
function F = FreeEnergy_var (ns, mu, V, mu0, V0, a, b, a0, b0, iFfx, iIn)

    F = 0;
    
    % Random effects
    % ---------------------------------------------------------------------
    % index of random effects
    iRfx = setdiff (1 : numel (mu), iFfx);
    iRfx = intersect (iRfx, iIn);
       
    % number of random effects
    n = length (iRfx);
    
    if n > 0
        % keep only random effects
        % '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        mu0 = mu0(iRfx);
        mu = mu(iRfx);
        V = V(iRfx, iRfx);
        V0 = V0(iRfx, iRfx);
        iv0 = VBA_inv (V0);
        a = a(iRfx);
        b = b(iRfx);
        a0 = a0(iRfx);
        b0 = b0(iRfx);

        % intermediary variables
        % '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        % expected group level variance
        E_lambda = a ./ b;
        % expected log group level variance
        E_log_lambda = psi (a) - log (b);
        % excursion
        e = mu - mu0;

        % compute free energy correction term
        % '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        % random effects
        F = F ...
            - 0.5 * ns * sum (log (E_lambda)) ...
            + sum ((a0 + 0.5 * ns - 1) .* E_log_lambda) ...
            - sum ((0.5 * ns * diag (V) + b0) .* E_lambda) ...
            + sum (a0 .* log (b0) + gammaln (b0)) ...
            - 0.5 * n * log (2 * pi) ...
            - 0.5 * VBA_logDet (V0) ...
            - 0.5 * e' * iv0 * e ...
            - 0.5 * trace (iv0 * V) ...
            + sum (VBA_entropy ('Gamma', a, 1 ./ b)) ...
            + VBA_entropy ('Gaussian', V);
    end
    
    % Fixed effects
    % ---------------------------------------------------------------------
    iFfx = intersect (iFfx, iIn);
    % number of random effects
    n = length (iFfx);
    if n > 0
        F = F + 0.5 * (ns - 1) .* n .* log (2 * pi);
    end
   
end

% check if parameter is fixed or random from group variance hyperparams
% =========================================================================
function il = isFixedEffect (a, b)
    il = isinf (a) & eq (b, 0);
end

% initialize within-subject priors
% =========================================================================
function priors = getSubjectPriors (priors_group, ind, nS)

    % start with fixed effects
    priors.muPhi = priors_group.muPhi;
    priors.SigmaPhi = nS * priors_group.SigmaPhi;
    
    priors.muTheta = priors_group.muTheta;
    priors.SigmaTheta = nS * priors_group.SigmaTheta;

    priors.muX0 = priors_group.muX0;
    priors.SigmaX0 = nS * priors_group.SigmaX0;
    
    % set random effects priors from group
    priors = updateSubjectPriors(priors, priors_group, ind);
end

% update within-subject priors from group posterior
% =========================================================================
% Update random effects only, as the prior for fixed effect is, by
% definition, fixed.
function priors = updateSubjectPriors (priors, posterior_grp, ind)
  
    priors.muPhi(ind.phi_rfx) = posterior_grp.muPhi(ind.phi_rfx);
    priors.SigmaPhi(ind.phi_rfx, ind.phi_rfx) = ...
        diag (posterior_grp.b_vPhi(ind.phi_rfx) ./ posterior_grp.a_vPhi(ind.phi_rfx));

    priors.muTheta(ind.theta_rfx) = posterior_grp.muTheta(ind.theta_rfx);
    priors.SigmaTheta(ind.theta_rfx, ind.theta_rfx) = ...
        diag (posterior_grp.b_vTheta(ind.theta_rfx) ./ posterior_grp.a_vTheta(ind.theta_rfx));

    priors.muX0(ind.x0_rfx) = posterior_grp.muX0(ind.x0_rfx);
    priors.SigmaX0(ind.x0_rfx, ind.x0_rfx) = ...
        diag (posterior_grp.b_vX0(ind.x0_rfx) ./ posterior_grp.a_vX0(ind.x0_rfx));
end

