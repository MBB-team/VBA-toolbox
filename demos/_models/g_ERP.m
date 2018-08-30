function [gx,dG_dX,dG_dPhi,d2G_dXdPhi] = g_ERP(x,Phi,u,inG)
% neural mass observation function (DCM for ERPs)

n = size(x,1);          % should be 9
nPhi = length(Phi);     % should be 1


g = exp(Phi(1)).*[0  0   0   0   0   0   0   0   1
                  0  0   0   0   0   0   0   0   0.3
                  0  0   0   0   0   0   0   0   -0.2
                  0  0   0   0   0   0   0   0   -0.7];


% state observation
gx = g*x;



%------ gradients evaluations ------%

% wrt the hidden states
dG_dX = g';


% wrt the observation parameters
dG_dPhi = (g*x)';


% mixed partial derivatives
d2G_dXdPhi(:,:,1) = zeros(n,nPhi);
d2G_dXdPhi(9,1,1) = exp(Phi(1));
d2G_dXdPhi(:,:,2) = 0.3*d2G_dXdPhi(:,:,1);
d2G_dXdPhi(:,:,3) = -0.2*d2G_dXdPhi(:,:,1);
d2G_dXdPhi(:,:,4) = -0.7*d2G_dXdPhi(:,:,1);
