function [priors,params2update] = VBA_fillInPriors(priors,dim,verbose)
% fills in priors structure with default priors if necessary
% function [priors,params2update] = VBA_fillInPriors(priors,dim,verbose)
% IN: [see VBA_check.m]
%   - priors: [can be left empty]
%   - dim: [required]
%   - verbose: flag for verbose mode
% OUT:
%   - priors: complete priors structure
%   - params2update: for VBA inversion

params2update.phi = 1:dim.n_phi;
params2update.theta = 1:dim.n_theta;
params2update.x0 = 1:dim.n;
for t=1:dim.n_t
    params2update.x{t} = 1:dim.n;
end
if ~isempty(priors)
    % Get user-specified priors
    fn = fieldnames(priors);
    priors0 = VBA_priors(dim,struct('binomial',0));
    fn0 = fieldnames(priors0);
    io = ismember(fn0,fn);
    ind = find(io==0);
    if ~isempty(ind)
        VBA_disp('Warning: could not find priors:',struct('verbose',verbose))
        for i = 1:length(ind)
            VBA_disp(['      - ',fn0{ind(i)}],struct('verbose',verbose));
            eval(['priors.',fn0{ind(i)},'=priors0.',fn0{ind(i)},';',])
        end
        VBA_disp('---> Using default (non-informative) priors',struct('verbose',verbose))
    end
    % check dimension and infinite precision priors
    if dim.n_theta > 0 % This finds which evolution params to update
        dpc = diag(priors.SigmaTheta);
        iz = find(dpc==0);
        if ~isempty(iz)
            params2update.theta = setdiff(1:dim.n_theta,iz);
        end
    end
    if dim.n_phi > 0 % This finds which observation params to update
        dpc = diag(priors.SigmaPhi);
        iz = find(dpc==0);
        if ~isempty(iz)
            params2update.phi = setdiff(1:dim.n_phi,iz);
        end
    end
    if dim.n > 0  % This finds which initial conditions to update
        dpc = diag(priors.SigmaX0);
        iz = find(dpc==0);
        if ~isempty(iz)
            params2update.x0 = setdiff(1:dim.n,iz);
        end
        for t=1:dim.n_t
            dpc = diag(priors.iQx{t});
            iz = find(isinf(dpc));
            if ~isempty(iz)
                params2update.x{t} = setdiff(1:dim.n,iz);
            end
        end
    end
    % insure vertical priors
    priors.muPhi = priors.muPhi(:);
    priors.muTheta = priors.muTheta(:);
    priors.muX0 = priors.muX0(:);
    if isfield(options.priors,'a_sigma')
        priors.a_sigma = options.priors.a_sigma(:);
        priors.b_sigma = options.priors.b_sigma(:);
    end
    
else % Build default (non-informative) priors
    priors = VBA_priors(dim,struct('binomial',0));
end

