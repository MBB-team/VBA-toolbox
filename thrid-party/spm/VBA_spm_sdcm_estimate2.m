function [DCM] = VBA_spm_sdcm_estimate2(P)
% Estimate parameters of a stochastic DCM for fMRI data (deconv HRF)
% FORMAT [DCM] = spm_sdcm_estimate2(DCM)   
%
% DCM  - the DCM or its filename
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau



% load DCM structure
if ~nargin
    P = VBA_spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    if isempty(P)
        return
    end
end
if isstruct(P)
    DCM = P;
    P   = ['DCM-' date];
else
    load(P)
end


disp(' ')
disp('-----------------------------')
disp('Inverting stochastic DCM...')


% Get data and remove confounds
y = DCM.Y.y';
X0 = DCM.Y.X0;
iX0 = pinv(X0'*X0)*X0';
nreg = size(y,1);
for i=1:nreg
    beta = iX0*y(i,:)';
    yc = X0*beta;
    y(i,:) = y(i,:) - yc';
end

% Unpack DCM specification
dt = DCM.U.dt;
uu = DCM.U.u';
[nu,nt] = size(uu);
[p,ny] = size(y);
A = DCM.a - eye(nreg);
for i=1:nu
    B{i} = DCM.b(:,:,i);
end
C = DCM.c;
if isfield(DCM,'d')
    for i=1:nreg
        try
            D{i} = DCM.d(:,:,i);
        catch
            D{i} = zeros(nreg,nreg);
        end
    end
else
    for i=1:nreg
        D{i} = zeros(nreg,nreg);
    end
end

% Prepare optional input for inversion routine
TR = DCM.Y.dt;
f_fname = @f_dcm4fmri;
g_fname = @g_Id;
microDT = 1e-1;             % micro-time resolution (in sec)
[inF] = prepare_dcm(A,B,C,D);
options.decim = max([1,floor(TR./microDT)]);
inF.deltat = TR./options.decim;
options.inF = inF;
options.microU = 1;
dim.n_theta = inF.indself;
dim.n_phi = 0;
dim.n = nreg;
dim.p = p;

% Build priors
alpha = 1e2/TR;
[priors] = getPriors(nreg,dim,alpha,ny);
options = getOptions(options,priors);

% Resample inputs on microtime integration grid
[u] = resampleU(uu,dt,nt,nu,ny,options);
u = [zeros(nu,1),u]; % add initial conditions

% Deconvolve fMRI data from HRF
[hrf] = VBA_spm_hrf(TR);
y0 = y;
y = [];
SNR = VBA_spm_getSNR(y0,uu,hrf);
for i=1:nreg
    try
        D = round(DCM.delays(i)./TR);
        y0i = [zeros(1,D),y0(i,1:end-D)];
    catch
        y0i = y0(i,:);
    end
    y(i,:) = VBA_spm_deconv(y0i,hrf,SNR);
end

% invert stochastic DCM using discrete-time VB routine:
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% extract relevant info and fills in DCM structure
[DCM] = exportDCMfromVBNLSS(posterior,out,DCM);

% change DCM output (account for HRF convolution)
for i=1:nreg
    y0i = conv(out.suffStat.gx(i,:),hrf);
    DCM.y(:,i) = y0i(1:end-length(hrf)+1)'; % predicted data
    DCM.R(:,i) = y0(i,:)' - DCM.y(:,i); % residuals of the model
end

disp('Inverting stochastic DCM... OK')
disp('-----------------------------')
disp(' ')

return



function [u] = resampleU(uu,dt,nt,nu,ny,options)
microDT = options.inF.deltat;
grid1 = 0:dt:dt*(nt-1);
grid2 = 0:microDT:microDT*options.decim*ny;
u = zeros(nu,ny);
for i=1:length(grid2)-1
    [tmp,ind1] = min(abs(grid2(i)-grid1));
    [tmp,ind2] = min(abs(grid2(i+1)-grid1));
    u(:,i) = mean(uu(:,ind1:ind2),2);
end
% [u,alpha] = spm_resample(full(uu),dt/microDT);



function [priors] = getPriors(nreg,dim,alpha,n_t)
priors.muX0 = 0.*ones(nreg,1);
priors.SigmaX0 = 0e-3*eye(nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.SigmaTheta = 1e-2*eye(dim.n_theta);
priors.muPhi = [];
priors.SigmaPhi = [];
% priors.SigmaTheta(options.inF.indself,options.inF.indself) = 0;

% NB on hyperpriors:
%   - fix state noise precision using high scale param
%   - use non-informative priors on the residual precision, with high
%   expectation.
% This is because of the first iteration of the hidden states posterior
% update, which has to deviate from the its prior predictive density (as
% derived from the deterministic inversion).
% The following iterations will then work with a realistic (expected)
% residual precision, and adapt.
SC = 1e0;
priors.a_alpha = SC*alpha;
priors.b_alpha = SC;
priors.a_sigma = 1e0;
priors.b_sigma = 1e-4;
for t = 1:n_t
    priors.iQx{t} = eye(dim.n);
    priors.iQy{t} = eye(dim.p);
end

function options = getOptions(options,priors)
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.gradF = 0;
options.updateHP = 1;
options.backwardLag = 2;
options.Laplace = 0;
% options.noSXi = 1;
% options.init0 = 0;
% options.embed = 0;


