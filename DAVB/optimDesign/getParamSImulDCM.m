function [theta,phi] = getParamSImulDCM(homogeneous,dim,options,pa,pc)
% builds the options structure from basic DCM info
%
% function [options,dim] = getOptions4dcm(A,B,C,D,TR,microDT,n_t,homogeneous)
% IN:
%   - homogeneous: flag for indicating whether the observation parameters
%   of the Ballon model are identical across ROIs
%   - dim: dimensions of the system
%   - options: the optional structure for VB inversion of the DCM
%   - pa/pc: values for the A and C simulated parameters
% OUT:
%   - options: optional structure for VB inversion of the specified model.

[nreg,nu] = size(options.inF.C);

% simulation  parameters
thetaHRF = [...
    0.4886
    0.6015
   -1.0861
   -0.8587
         0
    0.0545
    ];
phiHRF = [...
    2.1939
   -0.9883
   ];
thetaHRF = 0.*thetaHRF;
phiHRF = 0.*phiHRF;
% A matrix
na = length(find(options.inF.A(:)==1));
if numel(pa)==1
    t_A = repmat(pa,na,1);
else
    t_A = pa;
end
% self-inhibition gain
t_Aself = -0;
% B matrices
t_B{1} = [];
t_B{2} = [];
% C matrix
nc = length(find(options.inF.C(:)==1));
if numel(pa)==1
    t_C = repmat(pc,nc,1);
else
    t_C = pc;
end
% D matrices
t_D{1} = [];
t_D{2} = [];
t_D{3} = [];

%--- simulated evolution parameters: hemodynamic level
t_E0 = thetaHRF(1)*ones(nreg,1);       % HbO2 extraction fraction gain
t_tau0 = thetaHRF(2)*ones(nreg,1);     % mean blood transit time gain
t_kaf = thetaHRF(3)*ones(nreg,1);      % vasodilatory signal feedback regulation
t_kas = thetaHRF(4)*ones(nreg,1);      % vasodilatory signal decay gain
t_alpha = thetaHRF(6)*ones(nreg,1);    % vessel stifness gain

%--- simulated observation parameters
if ~homogeneous
    p_E0 = phiHRF(1)*ones(nreg,1);       % HbO2 extraction fraction gain
    p_epsilon = phiHRF(2)*ones(nreg,1);  % ratio of intra- and extravascular signal
else
    p_E0 = phiHRF(1);
    p_epsilon = phiHRF(2);
end

%--- Recollect paramters for simulated data
theta = zeros(dim.n_theta,1);
theta(options.inF.indA) = t_A;
for i=1:nu
    theta(options.inF.indB{i}) = t_B{i};
end
theta(options.inF.indC) = t_C;
for i=1:nreg
    theta(options.inF.indD{i}) = t_D{i};
end
theta(options.inF.indself) = t_Aself;
theta(options.inF.ind1) = t_E0;
theta(options.inF.ind2) = t_tau0;
theta(options.inF.ind3) = t_kaf;
theta(options.inF.ind4) = t_kas;
theta(options.inF.ind5) = t_alpha;
phi = zeros(dim.n_phi,1);
phi(options.inG.ind1) = p_E0;
phi(options.inG.ind2) = p_epsilon;


