function [fx,J,dfdp,d2fdxdp] = f_fullDCM4fmri(Xt,Theta,ut,inF)
% DCM for fMRI evolution function (attn: embedding of HRF!)
% function [fx,dF_dX,dF_dTheta] = f_fullDCM4fmri(Xt,Theta,ut,inF)
% This function evaluates the evolution funciton of the neuronal level in
% DCM for fMRI. For the sake of HRF convolution purposes, it also
% internally calls g_fullDCM4fmri.m so that the hemodynamic states are
% updated properly.


%- Neuronal states evolution
[fx,J,dfdp] = f_dcm4fmri(Xt,Theta,ut,inF);
d2fdxdp = [];

%- HRF convolution
g_fullDCM4fmri(Xt,[],[],inF.inG);
% This call is used to update the hemodynamic states that are convolved
% outputs of the neuronal states (Xt), as well as the necessary gradient,
% i.e. Jacobian and derivative w.r.t. hemodynamic parameters.


