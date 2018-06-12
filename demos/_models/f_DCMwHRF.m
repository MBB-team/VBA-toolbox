function [fx,J,dfdp] = f_DCMwHRF(Xt,Theta,ut,inF)
% DCM for fMRI evolution function (including Balloon model)
% function [fx,dF_dX,dF_dTheta] = f_DCMwHRF(Xt,Theta,ut,inF)
% This function evaluates the evolution function DCM for fMRI, including
% the Balloon HRF model.

ut = ut(inF.confounds.indu);

%- hidden states evolution
xn = Xt(inF.n5);
[fx,J,dfdp] = f_HRF3(Xt,Theta,xn,inF);
[fxn,Jn,dfndpn] = f_dcm4fmri(xn,Theta,ut,inF);

%- Reshape flow field and gradients
fx(inF.n5) = fxn;
J(inF.n5,inF.n5) = Jn;
if ~isempty(dfdp)
    dfdp(:,inF.n5)  = dfndpn;
end


