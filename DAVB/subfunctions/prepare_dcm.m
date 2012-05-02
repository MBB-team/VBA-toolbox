function [inF] = prepare_dcm(A,B,C,D)
% prepares matrices and indices vectors for DCM (without HRF)
% function [options] = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous)
% IN:
%   - A: binary matrix indicating where the connections are
%   - B: cell-array of binary matrices of modulatory effects
%   - C: binary matrix of input-state coupling
%   - D: cell-array of binay matrices for gating effects
% OUT:
%   - inF: the optional input structure to @f_dcm4fmri

[indA,indB,indC,indD,indself] = find_dcm(A,B,C,D);

inF.indself = indself;
inF.A = A;
inF.indA = indA;
inF.B = B;
inF.indB = indB;
inF.C = C;
inF.indC = indC;
inF.D = D;
inF.indD = indD;

[inF.dA,inF.dB,inF.dC,inF.dD] = get_dMatdvec(inF);


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

