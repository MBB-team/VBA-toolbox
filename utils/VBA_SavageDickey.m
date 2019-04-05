function [F2,po2] = VBA_SavageDickey(po1,pr1,F1,dim1,pr2)
% computes the free energy and posterior moments of a reduced model
% function [F2,po2] = VBA_SavageDickey(po1,pr1,F1,dim1,pr2)
% This function computes the approximate log model evidence and posterior
% density of a 'reduced' model, which is defined with respect to a 'full'
% model. This reduction is such that some of the non-zero prior variances
% of the full model are now set to zero. This is useful for testing
% assumptions of the form H0: theta=0. See Friston and Penny, 2011.
% IN:
%   - po1: posterior structure of the full model. This is typically
%   obtained from the inversion of the full model using
%   VBA_NLStateSpaceModel.m
%   - pr1: prior structure of the full model.
%   - F1: free energy of the full model.
%   - dim1: dimensions of the full model. Note: the dimensions of the
%   reduced model are identical, but some parameters along these dimensions
%   are now fixed to their prior mean.
%   - pr2: prior structure of the reduced model.
% OUT:
%   - F2: free energy of the reduced model.
%   - po2: posterior structure of the reduced model. note that the
%   posterior structure is obtained without having to invert the reduced
%   model...


pr1 = checkPriors(pr1,dim1);
pr2 = checkPriors(pr2,dim1);

F2 = F1;
if dim1.n_theta > 0
    mf = po1.muTheta;
    Sf = po1.SigmaTheta;
    mf0 = pr1.muTheta;
    Sf0 = pr1.SigmaTheta;
    mr0 = pr2.muTheta;
    Sr0 = pr2.SigmaTheta;
    [dF,mr,Sr] = VBA_spm_log_evidence(mf,Sf,mf0,Sf0,mr0,Sr0);
    po2.muTheta = mr;
    po2.SigmaTheta = Sr;
    F2 = F2 +dF;
end
if dim1.n_phi > 0
    mf = po1.muPhi;
    Sf = po1.SigmaPhi;
    mf0 = pr1.muPhi;
    Sf0 = pr1.SigmaPhi;
    mr0 = pr2.muPhi;
    Sr0 = pr2.SigmaPhi;
    [dF,mr,Sr] = VBA_spm_log_evidence(mf,Sf,mf0,Sf0,mr0,Sr0);
    po2.muPhi = mr;
    po2.SigmaPhi = Sr;
    F2 = F2 +dF;
end
if dim1.n > 0
    mf = po1.muX0;
    Sf = po1.SigmaX0;
    mf0 = pr1.muX0;
    Sf0 = pr1.SigmaX0;
    mr0 = pr2.muX0;
    Sr0 = pr2.SigmaX0;
    [dF,mr,Sr] = VBA_spm_log_evidence(mf,Sf,mf0,Sf0,mr0,Sr0);
    po2.muX0 = mr;
    po2.SigmaX0 = Sr;
    F2 = F2 +dF;
end

if isfield(po1,'a_sigma') && ~isempty(po1.a_sigma)
    mf = po1.a_sigma./po1.b_sigma;
    Sf = po1.a_sigma./(po1.b_sigma.^2);
    mf0 = pr1.a_sigma./pr1.b_sigma;
    Sf0 = pr1.a_sigma./(pr1.b_sigma.^2);
    mr0 = pr2.a_sigma./pr2.b_sigma;
    Sr0 = pr2.a_sigma./(pr2.b_sigma.^2);
    for iSource = 1:numel(po1.a_sigma)
        [dF,mr,Sr] = VBA_spm_log_evidence(mf(iSource),Sf(iSource),mf0(iSource),Sf0(iSource),mr0(iSource),Sr0(iSource));
        po2.b_sigma(iSource) = mr/Sr;
        po2.a_sigma(iSource) = po2.b_sigma(iSource)*mr;
        F2 = F2 +dF;
    end
end

if isfield(po1,'a_alpha') && ~isempty(po1.a_alpha) && ~isinf(po1.a_alpha)
    mf = po1.a_alpha./po1.b_alpha;
    Sf = po1.a_alpha./(po1.b_alpha^2);
    mf0 = pr1.a_alpha./pr1.b_alpha;
    Sf0 = pr1.a_alpha./(pr1.b_alpha^2);
    mr0 = pr2.a_alpha./pr2.b_alpha;
    Sr0 = pr2.a_alpha./(pr2.b_alpha^2);  
    [dF,mr,Sr] = VBA_spm_log_evidence(mf,Sf,mf0,Sf0,mr0,Sr0);
    po2.b_alpha = mr./Sr;
    po2.a_alpha = po2.b_alpha.*mr;
    F2 = F2 +dF;

end   

function priors = checkPriors(priors,dim)

options.sources = struct('type', 0 , ...
                         'out' , 1:dim.p          );
                     
fn                  = fieldnames(priors);
priors0             = VBA_defaultPriors(dim,options);
fn0                 = fieldnames(priors0);
io                  = ismember(fn0,fn);
ind                 = find(io==0);
if ~isempty(ind)
    for i = 1:length(ind)
        eval(['priors.',fn0{ind(i)},'=priors0.',fn0{ind(i)},';',])
    end
end




