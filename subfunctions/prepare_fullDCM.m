function [options] = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD,sources)
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

if nargin < 7
    homogeneous = 0;
else
    homogeneous = ~~homogeneous;
end

if ~exist('sources','var')
    sources(1).out=size(A,1);
    if exist('hA','var')
        sources(2).out=size(hA,1);
    end
end
extended = numel(sources)>1;

%- prepare neural evolution function parameters indices and matrices
[inF] = prepare_dcm(A,B,C,D);

dim.n = size(A,1);
dim.n_u = size(B,2);


%- prepare decoding function parameters indices and matrices
if extended
    [inF,nresp] = extend_dcm(inF,hA,hB,hC,hD,dim,sources);  
else
    nresp = 0;
end

%- define hemodynamic parameters indices
nreg=dim.n;
if extended
    offset = inF.indhself(end);
else
    offset = inF.indself;
end
inF.ind1 = [1:5:5*nreg] + offset;
inF.ind2 = [2:5:5*nreg] + offset;
inF.ind3 = [3:5:5*nreg] + offset;
inF.ind4 = [4:5:5*nreg] + offset;
inF.ind5 = [5:5:5*nreg] + offset;


%- define hidden states indices
inF.n1 = 1:5:5*nreg;
inF.n2 = 2:5:5*nreg;
inF.n3 = 3:5:5*nreg;
inF.n4 = 4:5:5*nreg;
inF.n5 = 5:5:5*nreg;

if extended
    try
        offset = inF.n5(end);
    catch
        offset=0;
    end
    inF.r =  (1:nresp) + offset;
end
inG.n1 = inF.n1;
inG.n2 = inF.n2;
inG.n3 = inF.n3;
inG.n4 = inF.n4;
inG.n5 = inF.n5;

if extended
    inG.r = offset+(1:nresp);        
    for i=2:numel(sources)
        sourceRespIdx{i-1} = sources(i).out - sources(1).out(end);
    end
    inG.sourceRespIdx = sourceRespIdx;
end


%- prepare observation function parameters indices and matrices
if ~homogeneous
    inG.ind1 = 1:2:2*nreg;
    inG.ind2 = 2:2:2*nreg;
else
    inG.ind1 = 1;
    inG.ind2 = 2;
end
%- define decoding parameters indices
if extended
    % default is one parameter for each source, which is a constant offset
    inG.indr = inG.ind2(end) + (1:numel(sources)-1);
end

% dimensions 
dim.n_theta = inF.ind5(end);
if extended
    dim.n_phi = inG.indr(end);
else
    dim.n_phi = inG.ind2(end);
end
dim.p = nreg + nresp ;
dim.n = 5*nreg+nresp;

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

inF.extended=extended;
inG.extended=extended;
options.extended=extended;


if extended
    options.sources=sources;
end

options.inF = inF;
options.inG = inG;

options.dim = dim;



