function [u,options,dim] = addConfounds2dcm(X0,u,options,dim)
% this function modifies the i/o of DCM inversion to cope with confounds
% function [u,options,dim] = addConfounds2dcm(X0,u,options,dim)
% This function adds in a dummy input that indexes time, which is used in
% the observation function of DCM to extract the appropriate line of the
% confounds matrix. Required additional parameters and indexes are passed
% to the evolution/observation functions through the inF/inG structures.
% IN:
%   - X0: n_tXn_p0 matrix of confounds (n_t is the nmber of time samples
%   in the data)
%   - u/options/dim: usual input to VBA inversion of a DCM.
% OUT:
%   - u/options/dim: input to VBA inversion of a DCM (augmented with
%   confounds).

if ~isempty(X0)
    nreg = size(options.inF.A,1);
    nu = size(u,1);
    [n_t,n_p0] = size(X0);
    
    S0 = options.priors.SigmaPhi;
    options.priors.muPhi = [options.priors.muPhi;zeros(nreg*n_p0,1)];
    options.priors.SigmaPhi = eye(dim.n_phi+nreg*n_p0);
    options.priors.SigmaPhi(1:dim.n_phi,1:dim.n_phi) = S0;
    options.inF.confounds.indu = 1:nu;
    options.inG.confounds.indu = 1:nu;
    options.inG.confounds.indt = nu+1;
    options.inG.confounds.X0 = X0;
    options.inG.confounds.indp = dim.n_phi+1:dim.n_phi+nreg*n_p0;
    ux0 = 1:n_t;
    if options.microU && ~isequal(options.decim,1)
        ux0 = kron(ux0,ones(1,options.decim));
    end
    u = [u;[ux0,zeros(size(u,2)-size(ux0,2))]];
    dim.n_phi = dim.n_phi + nreg*n_p0;
end

