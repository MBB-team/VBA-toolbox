function [options] = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous)
% precalculates intermediary variables for the VB inversion of DCM for fMRI
% function [options] = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous)
% IN:
%   - A: binary matrix indicating where the connections are
%   - B: cell-array of binary matrices of modulatory effects
%   - C: binary matrix of input-state coupling
%   - D: cell-array of binay matrices for gating effects
%   - TR: fMRI reptition time (data sampling period)
%   - microDT: micro-time resolution for ODE integration
%   - homogeneous: flag for indicating whether the observation parameters
%   of the Ballon model are identical across ROIs
% OUT:
%   - options: incomplete optinal structure for VB inversion of the
%   specified model (this does not include priors)...
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

if nargin < 7
    homogeneous = 0;
else
    homogeneous = ~~homogeneous;
end

%- prepare neural evolution function parameters indices and matrices
[inF] = prepare_dcm(A,B,C,D);
nreg = size(A,1);

%- define hemodynamic parameters indices
inF.ind1 = [1:5:5*nreg] + inF.indself;
inF.ind2 = [2:5:5*nreg] + inF.indself;
inF.ind3 = [3:5:5*nreg] + inF.indself;
inF.ind4 = [4:5:5*nreg] + inF.indself;
inF.ind5 = [5:5:5*nreg] + inF.indself;

%- define hidden states indices
nreg = size(A,1);
inF.n1 = 1:5:5*nreg;
inF.n2 = 2:5:5*nreg;
inF.n3 = 3:5:5*nreg;
inF.n4 = 4:5:5*nreg;
inF.n5 = 5:5:5*nreg;
inG.n1 = inF.n1;
inG.n2 = inF.n2;
inG.n3 = inF.n3;
inG.n4 = inF.n4;
inG.n5 = inF.n5;


%- prepare observation function parameters indices and matrices
if ~homogeneous
    inG.ind1 = 1:2:2*nreg;
    inG.ind2 = 2:2:2*nreg;
else
    inG.ind1 = 1;
    inG.ind2 = 2;
end

%- finalize options structure
options.decim = max([1,ceil(TR./microDT)]);
inF.deltat = TR./options.decim;
inF = orderfields(inF);
inF.confounds.indu = 1:size(C,2);
inF.confounds.indp = [];
inF.fullDCM = 1;
inF.linearized = 0;
inF.xshift = 0;
inF.logx2 = 1;
inG.fullDCM = 1;
inG.homogeneous = homogeneous;
inG.TE = 0.04;
inG.confounds.indu = 1:size(C,2);
inG.confounds.indt = [];
inG.confounds.X0 = [];
inG.confounds.indp = [];

options.inF = inF;
options.inG = inG;

