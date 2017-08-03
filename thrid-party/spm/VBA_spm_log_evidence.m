function [F,sE,sC] = VBA_spm_log_evidence(varargin)
% Return the log-evidence of a reduced model (under Laplace approximation)
% FORMAT [F,sE,sC] = spm_log_evidence(qE,qC,pE,pC,rE,rC)
% FORMAT [F,sE,sC] = spm_log_evidence(qE,qC,pE,pC,priorfun,varargin)
% FORMAT [F,sE,sC] = spm_log_evidence(qE,qC,pE,pC)
%
% qE,qC    - posterior expectation and covariance of full model
% pE,pC    - prior expectation and covariance of full model
% rE,rC    - prior expectation and covariance of reduced model
% or 
% priorfun - inline function that returns prior moments
%            {rE rC} = priorfun(varargin{:})
%
% or (if omitted) rE = 0 and rC = 0;
%
% F        - reduced log-evidence: ln p(y|reduced model) - ln p(y|full model)
% [sE,sC]  - posterior expectation and covariance of reduced model
%
%--------------------------------------------------------------------------
% This routine assumes the reduced model is nested within a full model and
% that the posteriors (and priors) are Gaussian. Nested here means that the
% prior precision of the reduced model, minus the prior precision of the
% full model is positive definite. We additionally assume that the prior
% means are unchanged. The two input argument formats are for use with
% spm_argmax.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_log_evidence.m 4281 2011-03-31 19:49:57Z karl $
 
% Compute reduced log-evidence
%==========================================================================
 
% check to see if priors are specified by a function
%--------------------------------------------------------------------------
qE = varargin{1};
qC = varargin{2};
pE = varargin{3};
pC = varargin{4};
try
    priors = varargin{5}(varargin{6:end});
    rE     = priors{1};
    rC     = priors{2};
catch
    try
        rE = varargin{5};
        rC = varargin{6};
    catch
        n  = size(qC,1);
        rE = sparse(n,1);
        rC = sparse(n,n);
    end
end
 
% reduced subspace 
%--------------------------------------------------------------------------
qE  = VBA_spm_vec(qE);
pE  = VBA_spm_vec(pE);
rE  = VBA_spm_vec(rE);
 
if nargout < 2
    dE  = pE - rE;
    dC  = pC - rC;
    k   = find(dE | any(dC,2));
    if ~isempty(k)
        qE  = qE(k);
        pE  = pE(k);
        rE  = rE(k);
        qC  = qC(k,k);
        pC  = pC(k,k);
        rC  = rC(k,k);
    else
        
        % the reduced and full models are the same
        %------------------------------------------------------------------
        F   = 0;
        return
    end
end

% fix tolerance for matrix inversions
%--------------------------------------------------------------------------
TOL   = exp(-16);

% remove fixed parameters under full model
%--------------------------------------------------------------------------
i     = find(diag(pC));

% preliminaries
%--------------------------------------------------------------------------
qP    = VBA_spm_inv(qC(i,i),TOL);
pP    = VBA_spm_inv(pC(i,i),TOL);
rP    = VBA_spm_inv(rC(i,i),TOL);
sP    = qP + rP - pP;
sC    = VBA_spm_inv(sP,TOL);
sE    = qP*qE(i) + rP*rE(i) - pP*pE(i);

% log-evidence
%--------------------------------------------------------------------------
F     = VBA_spm_logdet(rP*qP*sC*pC(i,i)) ...
      - (qE(i)'*qP*qE(i) + rE(i)'*rP*rE(i) - pE(i)'*pP*pE(i) - sE'*sC*sE);
F     = F/2;
    
% restore full conditional density
%--------------------------------------------------------------------------
if nargout > 1
    rE(i)   = sC*sE;
    rC(i,i) = sC;
    sE      = VBA_spm_unvec(rE,varargin{1});
    sC      = rC;
end
