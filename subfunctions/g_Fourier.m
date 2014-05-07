function [gx,dG_dX,dG_dPhi] = g_Fourier(Xt,Phi,ut,inG)
% Fourier basis set observation function (dummy HRF model)

% creating the discrete cosine Fourier set
K = size(Phi,1);     % number of Fourier basis functions
p = inG.p;           % number of time bins in the fMRi time series
X = zeros(p,K);
grid = [0:p-1]';
for k = 1:K
    X(:,k) = sin(grid.*(k-1)*pi./(p-1));
end
X(:,1) = 1;

% build data prediction (GLM)
gx = X*Phi;

% build the observation function gradients
dG_dX = zeros(1,length(gx));

dG_dPhi = X';






