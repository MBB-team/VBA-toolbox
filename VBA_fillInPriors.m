function [options,params2update] = VBA_fillInPriors(dim,options)
% fills in priors structure with default priors if necessary
% function [priors,params2update] = VBA_fillInPriors(priors,dim,verbose)
% IN: [see VBA_check.m]
%   - priors: [can be left empty]
%   - dim: [required]
%   - verbose: flag for verbose mode
% OUT:
%   - priors: complete priors structure
%   - params2update: for VBA inversion

% default indexing of unknown model variabes 
params2update.phi = 1:dim.n_phi;
params2update.theta = 1:dim.n_theta;
params2update.x0 = 1:dim.n;
for t=1:dim.n_t
    params2update.x{t} = 1:dim.n;
end

% VBA default priors
default_priors = VBA_priors(dim,options);
if ~isfield(options,'priors')
    options.priors = struct();
end

% fill in VBA's priors structure with default priors
options.priors = check_struct(options.priors,default_priors);



    
        


