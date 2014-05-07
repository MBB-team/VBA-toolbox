% Demo for DCM for fMRI (without Balloon observer)
% This demo inverts the DCM for fMRI model, without the ballon model, which
% is replaced by a nonlinear sigmoid observation function.


% close all
clear variables



%% First get neural level dynamics
n_t = 1.5e2;
TR = 2e0;
microDT = 1e-1;
nreg = 3;
nu = 2;
u       = zeros(2,n_t);         % inputs
u(1,2:15) = 1;
u(1,30:31) = 1;
u(1,50:51) = 1;
u(2,30:31) = 1;
u(2,60:61) = 1;
alpha   = Inf;
sigma   = Inf;
% DCM specification
A = [0 1 1
     1 0 1
     0 1 0];
B{1} = zeros(nreg,nreg);
B{2} = [0 0 0
        1 0 0
        0 0 0];
C = [1 0
     0 0
     0 0];
D{1} = [0 0 0
        0 0 0
        0 1 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);
% Prepare precalculated matrices for DCM inversion
[options.inF] = prepare_dcm(A,B,C,D);
inG.G0 = 50;
inG.k = 1;
% A matrix
t_A = exp([ -1.5
            -1.5
            -0.5
            -2.5
            -1.5 ]);
% self-inhibition gain
t_Aself = -1;
% B matrices
t_B{1} = [];
t_B{2} = exp([ -0.5 ]);
% C matrix
t_C = exp([ +0.1 ]);
% D matrices
t_D{1} = 0*exp([ -2 ]);
t_D{2} = [];
t_D{3} = [];
% fill in theta
theta(options.inF.indA) = t_A;
for i=1:nu
    theta(options.inF.indB{i}) = t_B{i};
end
theta(options.inF.indC) = t_C;
for i=1:nreg
    theta(options.inF.indD{i}) = t_D{i};
end
theta(options.inF.indself) = t_Aself;
theta = theta(:);
phi     = [];
% Build options structure for temporal integration of SDE
options.inG     = inG;
options.decim = max([1,floor(TR./microDT)]);
options.inF.deltat = TR./options.decim;
dim.n_theta         = length(theta);
dim.n_phi           = 0;
dim.n               = nreg;
dim.p               = nreg;
dim.n_t             = n_t;
f_fname = @f_dcm4fmri;
g_fname = @g_Id;
% get micro-time dynamics
out0.options = options;
out0.options.dim = dim;
out0.options.OnLine = 0;
out0.options.microU = 0;
out0.suffStat.dx = [];
out0.options.f_fname = f_fname;
out0.options.g_fname = g_fname;
[xx] = VBA_microTime(...
    struct('muTheta',theta,'muPhi',phi,'muX0',[0;0;0]),u,out0);


%% then get full dynamics (including hemodynamics)
[thetaHRF,phiHRF]   = get_HRFparams(TR,options.inF.deltat);
f_fname2             = @f_DCMwHRF;
g_fname2             = @g_HRF3;
[options]           = prepare_fullDCM(A,B,C,D,TR,options.inF.deltat);
options.inG.TE      = 0.04;
% simulated evolution parameters: hemodynamic level
t_E0 = thetaHRF(1)*ones(nreg,1);       % HbO2 extraction fraction gain
t_tau0 = thetaHRF(2)*ones(nreg,1);     % mean blood transit time gain
t_kaf = thetaHRF(3)*ones(nreg,1);      % vasodilatory signal feedback regulation
t_kas = thetaHRF(4)*ones(nreg,1);      % vasodilatory signal decay gain
t_alpha = thetaHRF(6)*ones(nreg,1);    % vessel stifness gain
% simulated observation parameters
p_E0 = phiHRF(1)*ones(nreg,1);       % HbO2 extraction fraction gain
p_epsilon = phiHRF(2)*ones(nreg,1);  % ratio of intra- and extravascular signal
% fill in theta
theta(options.inF.ind1) = t_E0;
theta(options.inF.ind2) = t_tau0;
theta(options.inF.ind3) = t_kaf;
theta(options.inF.ind4) = t_kas;
theta(options.inF.ind5) = t_alpha;
phi(options.inG.ind1) = p_E0;
phi(options.inG.ind2) = p_epsilon;
phi = phi(:);
dim.n_theta         = options.inF.ind5(end);
dim.n_phi           = 2*nreg;
dim.n               = 5*nreg;
dim.p               = nreg;
dim.n_t             = n_t;
X0 = kron(ones(nreg,1),[0;0;0;0;0]);
% get micro-time dynamics
out.options = options;
out.options.dim = dim;
out.options.OnLine = 0;
out.options.microU = 0;
out.suffStat.dx = [];
out.options.f_fname = f_fname2;
out.options.g_fname = g_fname2;
[xx0,y] = VBA_microTime(...
    struct('muTheta',theta,'muPhi',phi,'muX0',X0),u,out);


%% Finally, convolve neural dynamics with canonical HRF
f_fname3     = @f_HRF2;               % Ballon model evolution function
g_fname3     = @g_HRF3;               % Balloon model observation function
options.inF.fullDCM = 0;
options.inG.fullDCM = 0;
options.inG.TE      = 0.04;
options.microU      = 1;
dim.n_theta         = 6;
dim.n_phi           = 2;
dim.n               = 4;
dim.p               = 1;
dim.n_t             = n_t;
dim.u               = 1;
% get micro-time dynamics
out1.options = options;
out1.options.dim = dim;
out1.dim = dim;
out1.options.OnLine = 0;
out1.options.microU = 1;
out1.suffStat.dx = [];
out1.options.f_fname = f_fname3;
out1.options.g_fname = g_fname3;
yhrf = zeros(size(xx));
xx1 = zeros(size(xx0));
xx1(options.inF.n5,:) = xx;
for i=1:nreg
    u = [xx(i,:)];
    [xxhrf,yhrf(i,:)] = VBA_microTime(...
        struct('muTheta',thetaHRF,'muPhi',phiHRF,'muX0',[0;0;0;0]),u,out1);
    xx1(options.inF.n1(i),:) = xxhrf(1,:);
    xx1(options.inF.n2(i),:) = xxhrf(2,:);
    xx1(options.inF.n3(i),:) = xxhrf(3,:);
    xx1(options.inF.n4(i),:) = xxhrf(4,:);
end

