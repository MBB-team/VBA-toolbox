function [DCM] = VBA_spm_dcm_fmri_check(P)
% post-hoc diagnostics for DCM (bilinear or nonlinear) of fMRI data
% FORMAT [DCM] = spm_dcm_fmri_check(DM)
% DCM - DCM structure or its filename
%
% This routine provides some diagnostics to ensure model inversion has
% converged. It plots the predicted and observed responses over all regions
% and provides the coefficient of determination ? or percent variance
% explained. This should normally be above 10%. An abnormally low
% coefficient of determination is highlighted in red. Quantitatively, one
% would normally expect to see one or more extrinsic (between source)
% connections with the strength of 1/8 Hz or greater. If all the extrinsic
% posterior expectations are below this value, then this suggests a failure
% of convergence or that the data are very noisy (possibly due to using
% very small regions of interest to summarise regional responses). Finally,
% the posterior correlations among all parameters are shown in terms of a
% correlation matrix. The number of effective parameters estimated is
% reported in terms of the (KL) divergence between the posterior and
% prior densities over parameters. This is divided by the log of the
% number of observations, by appealing to the Bayesian information
% criterion. The divergence corresponds to complexity or Bayesian
% surprise. Normally, one would expect the posterior and prior to diverge
% in a non-trivial fashion.
%
% Posterior densities are shown as bars with 90% confidence intervals in
% pink. An informed model inversion would normally provide posterior
% densities with confidence intervals that are, for some connections,
% displaced from prior expectations (at or around zero).
%
% This routine is compatible with DCM8, DCM10 and DCM12 files.
%__________________________________________________________________________
% Copyright (C) 20012 Wellcome Trust Centre for Neuroimaging
% Karl Friston
 
% $Id: spm_dcm_fmri_check.m karl $
 
%-Load DCM structure
%--------------------------------------------------------------------------
if ~nargin
  uiopen('load');
elseif isstruct(P)
  DCM = P;
else
  load(P)
end
 
% Assemble diagnostics
%==========================================================================
% coefficient of determination (percent variance explained)
%--------------------------------------------------------------------------
PSS  = sum(sum(DCM.y.^2));
RSS  = sum(sum(DCM.R.^2));
D(1) = 100*PSS/(PSS + RSS);
 
% largest absolute posterior expectation (extrinsic connections)
%--------------------------------------------------------------------------
try
   A = DCM.Ep.A;
catch
   A = DCM.A;
end
D(2) = max(max(abs(A - diag(diag(A)))));
 
% complexity and effective number of parameters estimated
%--------------------------------------------------------------------------
qE = VBA_spm_vec(DCM.Ep);
pE = VBA_spm_vec(DCM.M.pE);
qC = DCM.Cp;
pC = DCM.M.pC;
k  = rank(full(pC));
pC = VBA_spm_inv(pC);
 
D(3) = (trace(pC*qC) + (pE - qE)'*pC*(pE - qE) - VBA_spm_logdet(qC*pC) - k)/2;
D(3) = D(3)/log(DCM.v);
 
% Plot summary of inversion
%==========================================================================
VBA_spm_figure('GetWin','DCM diagnostics'); clf
 
% plot predicted and observed regional responses
%--------------------------------------------------------------------------
subplot(2,1,1);
t = (1:DCM.v)*DCM.Y.dt;
 
plot(t,DCM.y,t,DCM.y + DCM.R,':');
str = sprintf('variance explained %0.0f%%', D(1));
str = {'responses and predictions',str};
if D(1) > 10
  title(str,'FontSize',16);
else
  title(str,'FontSize',16,'Color','r');
end
xlabel('time {seconds}');
 
% posterior densities over A parameters
%--------------------------------------------------------------------------
try
  i = VBA_spm_fieldindices(DCM.Ep,'A');
catch
  i = 1 + (1:DCM.n^2);
end
qE = VBA_spm_vec(DCM.Ep);
qC = DCM.Cp;
subplot(2,2,3)
VBA_spm_plot_ci(qE(i),qC(i,i)), hold on
str = sprintf('largest connection strength %0.2f', D(2));
str = {'intrinsic and extrinsic connections',str};
if D(2) > 1/8
  title(str,'FontSize',16);
else
  title(str,'FontSize',16,'Color','r');
end
xlabel('parameter}');
axis square
 
% posterior correlations among all parameters
%--------------------------------------------------------------------------
subplot(2,2,4)
imagesc(cov2corr(DCM.Cp))
title('posterior correlations','FontSize',16)
str = sprintf('estimable parameters %0.0f', D(3));
str = {'posterior correlations',str};
if D(3) > 1
  title(str,'FontSize',16);
else
  title(str,'FontSize',16,'Color','r');
end
axis square