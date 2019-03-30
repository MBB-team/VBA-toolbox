function [results]=VBA_Shapley(posterior,out,varargin)
% [results]=VBA_Shapley(posterior,out,[options])
%
% Compute the Shapley values of the model's factors (inputs or paramters). 
% These scores measure the relative influence of these factor on the variaince 
% explained by the model.
% -------------------------------------------------------------------------
% IN:
%  - posterior, out: model structures
%  - varargin:       list of option/values pairs (see below)
% OUT:
%  - results: structure of results with the following fields, depending on
%  the chose coalitions
%  .parameters/inputs:  
%         Shapley values (first order) of paramters/inputs. This is a n x p
%         matrix where n is the number of factors of interest and p the 
%         number of observations in the model.
%  - interaction: Interaction score, ie relative change in Shapley values for a
%         perturbation in respective inputs. It is an array with the same
%         size as the inputs whose elements are similar to sv.
% -------------------------------------------------------------------------
% options:
%   coalitions:         ('parameters') | 'inputs' | 'interactions'
%      > factor of interest of the analysis. For 'interaction', compute the
%      Shapley value of parameters and their relative change for selective
%      input perturbations.
%   inputPerturbation   -> ('zero') | 'average' | 'average_nonzero'
%      > type of perturbation to apply on inputs: set to zero, average
%      accross the experiment, or set non zero inputs to their average
%   paramType           -> ('phi')  | 'theta'
%      > parameters of interest. 
%   paramIdx            -> (all)  | paramIdx
%      > restrict parameters of interest to the index list paramIdx
%   inputIdx            -> (all)  | inputIdx
%      > restrict inputs of interest to the index list inputIdx
%   obsIdx              -> (all)  | ObsIdx
%      > restrict observations of interest to the index list obsIdx
% -------------------------------------------------------------------------

%% Complete option structure
% -------------------------------------------------------------------------

if numel(varargin) == 1 && isstruct(varargin{1})
    options = varargin{1};
else
    options.coalitions = {'parameters','inputs','interactions'};
    options.inputPerturbation = {'zero','average','average_nonzero'};
    options.paramType = {'phi','theta'};
    options.paramIdx = 1:out.dim.n_phi ;
    options.inputIdx = 1:out.dim.u ;
    if size(out.y,2) == 1 % catch vertical data as unique observation
        options.obsIdx = 0;
    else
        options.obsIdx = 1:out.dim.p ;
    end
    parser = inputParser;
    parser.parse (varargin{:});
    options = parser.Results;
    %options = parseargs(options,varargin{:});
end

%% Prepare perturbation scheme
% -------------------------------------------------------------------------

% dimensions
nu = numel(options.inputIdx);
nw = numel(options.paramIdx);
if options.obsIdx==0 % catch vertical data as unique observation
    nResps=1;
else
    nResps = numel(options.obsIdx);
end

% factorial perturbations on coalitions of interest
switch options.coalitions
    case 'interactions'
        % factorial perturbation of parmeeters with normal inputs,
        kw = full(VBA_spm_perm_mtx(nw));
        k = [kw ones(2^nw,nu)] ;
        % factorial perturbation of parameters with each input respectively
        % pertrubed 
        ku = ones(nu)-eye(nu);
        for i=1:nu
            k = vertcat(k, [kw repmat(ku(i,:),size(kw,1),1)]); 
        end
        % plus factorial perturbation of inputs alone
        % (ie 2^nw + nu x 2^nw + 2^nu coalitions)
        k = vertcat(k,[ones(2^nu,nw) full(VBA_spm_perm_mtx(nu))]);
    case 'parameters'
        % factorial perturbation of paramters, normal inputs
        k = full(VBA_spm_perm_mtx(nw));
        k = [k ones(size(k,1),nu)];
    case 'inputs'
        % factorial perturbation of inputs, normal paramters
        k = full(VBA_spm_perm_mtx(nu));
        k = [ones(size(k,1),nw) k];
end
nk = size(k,1);

%% Compute explained variances
% -------------------------------------------------------------------------

% loop over coalitions
ve = nan(nk,nResps);
parfor t = 1:size(k,1)
    kt = k(t,:);
    w_perm = kt(1:nw);
    u_perm = kt(nw+(1:nu));
    ve(t,:) = explainedVar(posterior,out,options,u_perm,w_perm) ;
end

% normalize
ve1 = ve(1,:);
ve0 = ve(end,:);
for i=1:nResps
    ve(:,i) = (ve(:,i) - ve0(i) )/(ve1(i)-ve0(i)) ;
end

%% Compute Shapley values
% -------------------------------------------------------------------------

% restrict coalitions to effects of interest
switch options.coalitions
    case 'interactions'
        spl = [2^nw*ones(nu+1,1); 2^nu];
        k = mat2cell(k,spl,nw+nu);
        for i=1:nu+1
            k{i} = k{i}(:,1:nw);
        end
        k{nu+2} = k{nu+2}(:,nw+(1:nu));
        ve = mat2cell(ve,spl,nResps);
    case 'parameters'
        k = {k(:,1:nw)};
        ve = {ve};
    case 'inputs'
        k = {k(:,nw+(1:nu))};
        ve = {ve};      
end

% compute first order scores
for ii = 1:numel(k)
    n = size(k{ii},2);
    nn = factorial(n);
%     v{ii} = nan(n,nResps);
    for m=1:n % loop over players
        % Shapley coeficients
        i = k{ii}(:,m);
        z = sum(k{ii},2);
        coef = (2*i-1).*factorial(z-i).*factorial(n-z-(1-i))/nn;
        % compute shapley values per se
        v{ii}(m,:) =  coef'*ve{ii};
    end
end

% compute interactions if necessary
sv = v{1};
if strcmp(options.coalitions,'interactions')
    for iu=1:nu
        svi{iu} = (v{1}-v{iu+1})./v{1};
    end
    % shapley value of inputs
    svu = v{nu+2};
else
    svi={};
end

%% Store results
% -------------------------------------------------------------------------
switch options.coalitions
    case 'parameters'
        results.parameters = sv;
    case 'inputs'
        results.inputs = sv;
    case 'interactions'
        results.parameters = sv;
        results.inputs = svu;
        results.interactions = svi;
end

        
end
%% Subfunctions
% =========================================================================

% -------------------------------------------------------------------------
% Compute explained variance of the model given by posterior and out, with
% the inputs and paramters pertrubed according to u_swicth and w_switch
% respectively (switch = 1 -> normal, switch=0 -> perturbation)
% -------------------------------------------------------------------------
function v=explainedVar(posterior,out,options,u_switch,w_switch)

    % prevent unecessary bells and whistles
    out.options.verbose = 0;
    out.options.DisplayWin = 0;
    out.options.inF{1}.fast = true;

    % prepare degraded model
    % .....................................................................   
    
    % == pertub inputs
    % index of inputs to perturb
    inputIdx = options.inputIdx(u_switch==0);
    % apply perturbation
    switch options.inputPerturbation
        case 'zero'
            out.u(inputIdx,:) = 0;
        case 'average'
            out.u(inputIdx,:) = mean(out.u(inputIdx,:),2);
        case 'average_nonzero'
            for iu=inputIdx
                idxNZ = find(out.u(iu,:)~=0); 
                out.u(iu,idxNZ) = mean(out.u(iu,idxNZ));
            end
    end

    % == perturb parameters
    paramIdx = options.paramIdx(w_switch==0);
    switch options.paramType
        case 'phi'
            posterior.muPhi(paramIdx) = 0*posterior.muPhi(paramIdx);
        case 'theta'
            posterior.muTheta(paramIdx) = 0*posterior.muTheta(paramIdx);
    end

    % predict data
    % .....................................................................   
    [yp,~,~,~,er] = VBA_simulate (...
        out.options.dim.n_t,...
        out.options.f_fname,...
        out.options.g_fname,...
        posterior.muTheta,...
        posterior.muPhi,...
        out.u,...
        Inf,...
        Inf,...
        out.options,...
        posterior.muX0);
    g = yp-er;
    y = out.y;
    
    % if vertical data, transpose everything
    if options.obsIdx == 0 
        g = g';
        y = y';
        options.obsIdx = 1;
        out.options.isYout = out.options.isYout';
    end

    % compute explained variance
    % .....................................................................     
    v = nan(1,numel(options.obsIdx));
    for i=1:numel(options.obsIdx) % for each observation of interest
        obsIdx = options.obsIdx(i);
        in_idx = find(out.options.isYout(obsIdx,:) == 0);
        v(i) = 1-((var(y(obsIdx,in_idx)-g(obsIdx,in_idx))/var(y(obsIdx,in_idx))))  ;
    end

end

