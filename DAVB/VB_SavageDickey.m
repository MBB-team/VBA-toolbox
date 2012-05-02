function [F2,po2] = VB_SavageDickey(po1,pr1,F1,dim1,pr2)
% computes the free energy and posterior moments of a reduced model
% function [F2,po2] = VB_SavageDickey(po1,pr1,F1,dim1,pr2)
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
    if isequal(pr1.muTheta,pr2.muTheta) ...
            && isequal(full(pr1.SigmaTheta),full(pr2.SigmaTheta))
        po2.muTheta = po1.muTheta;
        po2.SigmaTheta = po1.SigmaTheta;
    else
        dP = diag(pr2.SigmaTheta);
        iout = find(dP==0);
        Pi1 = VB_inv(pr1.SigmaTheta);
        Pi2 = VB_inv(pr2.SigmaTheta);
        P1 = VB_inv(po1.SigmaTheta);
        P2 = P1;
        P2(iout,:) = 0;
        P2(:,iout) = 0;
        po2.SigmaTheta = VB_inv(P2);
        po2.muTheta = po2.SigmaTheta*(...
            P1*po1.muTheta + Pi2*pr2.muTheta - Pi1*pr1.muTheta );
        F2 = F2 ...
            + 0.5*VBA_logDet(Pi2) ...
            + 0.5*VBA_logDet(P1) ...
            - 0.5*VBA_logDet(P2) ...
            - 0.5*VBA_logDet(Pi1) ...
            - 0.5*po1.muTheta'*P1*po1.muTheta ...
            - 0.5*pr2.muTheta'*Pi2*pr2.muTheta ...
            + 0.5*pr1.muTheta'*Pi1*pr1.muTheta ...
            + 0.5*po2.muTheta'*P2*po2.muTheta;
    end
end
if dim1.n_phi > 0
    if isequal(pr1.muPhi,pr2.muPhi) ...
            && isequal(full(pr1.SigmaPhi),full(pr2.SigmaPhi))
        po2.muPhi = po1.muPhi;
        po2.SigmaPhi = po1.SigmaPhi;
    else
        iout = find(diag(pr2.SigmaPhi)==0);
        Pi1 = VB_inv(pr1.SigmaPhi);
        Pi2 = VB_inv(pr2.SigmaPhi);
        P1 = VB_inv(po1.SigmaPhi);
        P2 = P1;
        P2(iout,:) = 0;
        P2(:,iout) = 0;
        po2.SigmaPhi = VB_inv(P2);
        
        po2.muPhi = po2.SigmaPhi*(...
            P1*po1.muPhi + Pi2*pr2.muPhi - Pi1*pr1.muPhi );
%         F20 = F2 ...
%             + 0.5*VBA_logDet(Pi2) ...
%             + 0.5*VBA_logDet(P1) ...
%             - 0.5*VBA_logDet(P2) ...
%             - 0.5*VBA_logDet(Pi1) ...
%             - 0.5*po1.muPhi'*P1*po1.muPhi ...
%             - 0.5*pr2.muPhi'*Pi2*pr2.muPhi ...
%             + 0.5*pr1.muPhi'*Pi1*pr1.muPhi ...
%             + 0.5*po2.muPhi'*P2*po2.muPhi;
        
        
        fixed = setdiff(iout,find(diag(po1.SigmaPhi)==0));
        dmu = po1.muPhi(fixed) - pr2.muPhi(fixed);
        dmu0 = pr1.muPhi(fixed) - pr2.muPhi(fixed);
        F2 = F2 ...
            - 0.5*VBA_logDet(po1.SigmaPhi(fixed,fixed)) ...
            - 0.5*dmu'*VB_inv(po1.SigmaPhi(fixed,fixed))*dmu ...
            + 0.5*VBA_logDet(pr1.SigmaPhi(fixed,fixed)) ...
            + 0.5*dmu0'*VB_inv(pr1.SigmaPhi(fixed,fixed))*dmu0;
%         F20-F2
%         pause
        
        
    end
end
if dim1.n > 0
    if isequal(pr1.muX0,pr2.muX0) ...
            && isequal(full(pr1.SigmaX0),full(pr2.SigmaX0))
        po2.muX0 = po1.muX0;
        po2.SigmaX0 = po1.SigmaX0;
    else
        dP = diag(pr2.SigmaX0);
        iout = find(dP==0);
        Pi1 = VB_inv(pr1.SigmaX0);
        Pi2 = VB_inv(pr2.SigmaX0);
        P1 = VB_inv(po1.SigmaX0);
        P2 = P1;
        P2(iout,:) = 0;
        P2(:,iout) = 0;
        po2.SigmaX0 = VB_inv(P2);
        po2.muX0 = po2.SigmaX0*(...
            P1*po1.muX0 + Pi2*pr2.muX0 - Pi1*pr1.muX0 );
        F2 = F2 ...
            + 0.5*VBA_logDet(Pi2) ...
            + 0.5*VBA_logDet(P1) ...
            - 0.5*VBA_logDet(P2) ...
            - 0.5*VBA_logDet(Pi1) ...
            - 0.5*po1.muX0'*P1*po1.muX0 ...
            - 0.5*pr2.muX0'*Pi2*pr2.muX0 ...
            + 0.5*pr1.muX0'*Pi1*pr1.muX0 ...
            + 0.5*po2.muX0'*P2*po2.muX0;
    end
end


function priors = checkPriors(priors,dim)
fn                  = fieldnames(priors);
priors0             = VBA_priors(dim,struct('binomial',0));
fn0                 = fieldnames(priors0);
io                  = ismember(fn0,fn);
ind                 = find(io==0);
if ~isempty(ind)
    for i = 1:length(ind)
        eval(['priors.',fn0{ind(i)},'=priors0.',fn0{ind(i)},';',])
    end
end


