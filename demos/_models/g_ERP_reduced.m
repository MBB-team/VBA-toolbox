function [gx,dG_dX,dG_dPhi,d2G_dXdPhi] = g_ERP_reduced(x,Phi,u,inG)
% reduced neural mass observation function (DCM for ERPs)

n = size(x,1);          % should be 2
nPhi = length(Phi);     % should be 1



g = exp(Phi(1))*eye(2);
% state observation
gx = g*x;



%------ gradients evaluations ------%

% wrt the hidden states
dG_dX = g';


% wrt the observation parameters
dG_dPhi = (g*x)';


% mixed partial derivatives
d2G_dXdPhi(:,:,1) = [exp(Phi(1));0];
d2G_dXdPhi(:,:,2) = [0;exp(Phi(1))];

