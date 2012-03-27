function [gx,dG_dX,dG_dPhi,d2G_dXdPhi] = g_Fourier(Xt,Phi,ut,inG)
% Fourier basis set observation function (dummy HRF model)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

% creating the discrete cosine Fourier set
K = size(Phi,1);     % number of Fourier basis functions
p = inG.p;           % number of time bins in the fMRi time series
X = zeros(p,K);
grid = [0:p-1]';
for k = 1:K
    X(:,k) = sin(grid.*(k-1)*pi./(p-1));
end
X(:,1) = 1;

% figure,plot(X)
% pause


% build data prediction (GLM)
gx = X*Phi;

% build the observation function gradients
dG_dX = zeros(1,length(gx));

dG_dPhi = X';

for n = 1:length(gx)
    d2G_dXdPhi(:,:,n) = zeros(1,length(Phi));
end





