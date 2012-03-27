function [options] = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous)
% precalculates intermediary variables for DCM for fMRI VB inversion

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

%- prepare evolution function parameters indices and matrices
inF.A = A;
inF.B = B;
inF.C = C;
inF.D = D;
[inF.indA,inF.indB,inF.indC,inF.indD,inF.indself] = find_dcm(A,B,C,D);
[inF.dA,inF.dB,inF.dC,inF.dD] = get_dMatdvec(inF);
inF.ind1 = [1:5:5*nreg] + inF.indself;
inF.ind2 = [2:5:5*nreg] + inF.indself;
inF.ind3 = [3:5:5*nreg] + inF.indself;
inF.ind4 = [4:5:5*nreg] + inF.indself;
inF.ind5 = [5:5:5*nreg] + inF.indself;

%- prepare observation function parameters indices and matrices
if ~homogeneous
    inG.ind1 = 1:2:2*nreg;
    inG.ind2 = 2:2:2*nreg;
else
    inG.ind1 = 1;
    inG.ind2 = 2;
end

%- build options structure
inF.fullDCM = 1;
inF.linearized = 0;
% inF.xshift = 1e-1;
inG.fullDCM = 1;
inG.homogeneous = homogeneous;
inG.TE = 0.04;
options.decim = max([1,ceil(TR./microDT)]);
inF.deltat = TR./options.decim;

inF = orderfields(inF);
options.inF = inF;
options.inG = inG;



function [dA,dB,dC,dD] = get_dMatdvec(inF)
if ~isempty(inF.indA)
    dA = dMatdvec(inF.A);
else
    dA = [];
end
dB = cell(length(inF.B),1);
for i=1:length(inF.B)
    if ~isempty(inF.indB{i})
        dB{i} = dMatdvec(inF.B{i});
    else
        dB{i} = [];
    end
end
if ~isempty(inF.indC)
    dC = dMatdvec(inF.C);
else
    dC = [];
end
dD = cell(length(inF.D),1);
for i=1:length(inF.D)
    if ~isempty(inF.indD{i})
        dD{i} = dMatdvec(inF.D{i});
    else
        dD{i} = [];
    end
end

function C = dMatdvec(A)
A = ~~A;
ind = find(A~=0);
n = numel(A);
ni = length(ind);
C = zeros(n,ni);
for i=1:length(ind)
    C(ind(i),i) = 1;
end

function [indA,indB,indC,indD,indself] = find_dcm(A,B,C,D)
ia = find(A~=0);
if ~isempty(ia)
    indA = 1:length(ia);
else
    indA = [];
end
ib = [];
indB = cell(length(B),1);
for i=1:length(B)
    tmp = find(B{i}~=0);
    if ~isempty(tmp)
        indB{i} = length(ia)+length(ib)+1:...
            length(ia)+length(ib)+length(tmp);
    else
        indB{i} = [];
    end
    ib = [ib;tmp];
end
ic = find(C~=0);
if ~isempty(ic)
    indC = length(ia)+length(ib)+1:length(ia)+length(ib)+length(ic);
else
    indC = [];
end
id = [];
indD = cell(length(D),1);
for i=1:length(D)
    tmp = find(D{i}~=0);
    if ~isempty(tmp)
        indD{i} = length(ia)+length(ib)+length(ic)+length(id)+1:...
            length(ia)+length(ib)+length(ic)+length(id)+length(tmp);
    else
        indD{i} = [];
    end
    id = [id;tmp];
end
indself = length(ia)+length(ib)+length(ic)+length(id)+1;
