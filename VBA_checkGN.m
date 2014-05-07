function [dx,flag] = VBA_checkGN(S,dx0)
% Levenberg-Marquardt regularization on Gauss-Newton update step
% [dx] = VBA_checkGN(S,dx0)
% This function checks the Gauss-Newton update required for the
% optimization of the variational energies of parameters. This
% regularization should not happen in normal cases (nonlinear Gaussian
% hierarchical models).
% IN:
%   - S: the nxn covariance matrix at the current mode
%   - dx0: the Gauss-Newton update step for the mode
% OUT:
%   - dx: the (Levenberg-Marquardt) regularized update step
%   - flag=1 if regularization, flag=0 if not.

n = numel(S);
% Compute smallest eigenvalue
if n == 1
    sev = S;
else
%     warning('off','all')
    sev = eigs(S,1,'sr',struct('disp',0));
%     warning('on','all')
end
if sev <= 0
    [V,D] = eig(S);
    D = -diag(D).^-1;
    tau = 1./max(abs(D));
    P = V*diag(exp(tau.*D))*V';
    dx = dx0 - P*dx0;
    flag = 1;
else
    dx = dx0;
    flag = 0;
end