function [DJS,b,my,Vy,mu,Q,pm] = dcm_risk(dcm_fnames,pm,split,fam,P)
% derive the Laplace-Chernoff risk for a given comparison set
% IN:
%   - dcm_fnames: nmX1 cell array of file names containing the DCM
%   structures
%   - pm: nmX1 vector of model prior probabilities {ones(nm,1)/nm}
%   - split: flag for split-Laplace derivation of the prior predictive
%   densities {0}
%   - fam: nfX1 cell array of model indices defining a partition of the
%   comparison set. This is used to derive the Laplace-Chernoff risk
%   of the associated family inference; When left empty (default),
%   inference is at the model level {[]}
%   - P: 2X1 vector prior means of evolution parameters (P(1)) and
%   precision hyperparameters (P(2)) {[1e-2;5e-2]}
% OUT:
%   - DJS: Jensen-Shannon divergence between the prior predictive densities
%   of the models/families included in the comparison set
%   - b: Laplace-Chernoff upper bound on the model selection error rate
%   - my: 1st-order moment of the marginal prior predictive density
%   - Vy: 2nd-order moment of the marginal prior predictive density
%   - mu: cell array of 1st-order moments of the conditional (upon
%   models/families) prior predictive densities
%   - Q: cell array of 2nd-order moments of the conditional (upon
%   models/families) prior predictive densities
%   - pm: vector of prior probabilities of models/families
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

nm = length(dcm_fnames);
try; pm=pm./sum(pm); catch; pm=ones(nm,1)./nm; end
try; split; catch; split = 0; end
try; fam; catch; fam = []; end
try; P; catch; P = [1e-2;5e-2]; end

if split
    disp('Deriving the split-Laplace approximation to the prior predictive densities:')
else
    disp('Deriving the Laplace approximation to the prior predictive densities:')
end

mu = cell(0);
Q = cell(0);
for i=1:nm
    fprintf(1,[' - model ',num2str(i),'/',num2str(nm),'...'])
    if isstruct(dcm_fnames{i})
        DCM = dcm_fnames{i};
    elseif exist(dcm_fnames{i},'file')
        load(dcm_fnames{i})
    else
        disp([' ! Warning: ',dcm_fnames{i},' does not exist!'])
        return
    end
    % export to VBA format
    [y,u,f_fname,g_fname,dim,options] = dcm2vba(DCM,0);
    % shift prior mean of connectivity parameters
    options.priors.muTheta(1:options.inF.indself) = P(1);
    % inflate noise variance
    options.priors.a_sigma = P(2);
    options.priors.b_sigma = 1e0;
    if split
        % get split-Laplace approx to the prior predictive density
        % !! do not split hemodynamic parameters !!
        nosplit.phi = [1:dim.n_phi]';
        nosplit.theta = [options.inF.indself:dim.n_theta]';
        nosplit.x0 = [1:dim.n]';
        [mu{end+1},Q{end+1}] = ...
            splitLaplace(u,f_fname,g_fname,dim,options,2,nosplit);
    else
        % get Laplace approx to the prior predictive density
        [mu{end+1},Q{end+1}] = ...
            getLaplace(u,f_fname,g_fname,dim,options,0);
    end
    fprintf(1,[' OK.'])
    fprintf(1,'\n')
    
end

% if family inference, pool prior predictive densities
if ~isempty(fam)
    fprintf(1,'Pooling prior predictive densities for family inference...')
    [mu,Q,pm] = pool_Laplace(mu,Q,fam);
    fprintf(1,[' OK.'])
    fprintf(1,'\n')
end

% derive approximation to the Chernoff bound
if length(mu) > 1
    fprintf(1,'Deriving the Laplace-Chernoff risk...')
    [DJS,b,my,Vy] = JensenShannon(mu,Q,pm,0,'2');
    fprintf(1,[' OK.'])
    fprintf(1,'\n')
else
    DJS = [];
    b = [];
    my = [];
    Vy = [];
end




