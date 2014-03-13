function [options] = prepare_fullDCM_old(A,B,C,D,deltat,n_t,decim,priors)
% prepares matrices and indices vectors for DCM for fMRI
clear g_fullDCM4fmri
clear f_fullDCM4fmri

%- Prepare observation function preprocessed inputs
nreg = size(A,1);
for i=1:nreg
    inG.indHRF{i} = 6*(i-1)+1:6*i;
end
inG.n_t = n_t;
inG.deltat = deltat;
inG.n = nreg;
inG.n1 = 1:4:4*nreg;
inG.n2 = 2:4:4*nreg;
inG.n3 = 3:4:4*nreg;
inG.n4 = 4:4:4*nreg;
inG.ind1 = 1:6:6*nreg;
inG.ind2 = 2:6:6*nreg;
inG.ind3 = 3:6:6*nreg;
inG.ind4 = 4:6:6*nreg;
inG.ind5 = 5:6:6*nreg;
inG.ind6 = 6:6:6*nreg;
inG.fullDCM = 1;
inG.Phi = priors.muPhi;

%- Prepare evolution function preprocessed inputs
[indA,indB,indC,indD,indself] = find_dcm(A,B,C,D);
inF.deltat = deltat;
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
inF.inG = inG;

% Build options structure
options.priors      = priors;
options.annealing   = 0;
options.ignoreMF    = 1;
options.inF     = inF;
options.inG     = inG;
% options.u0      = 0*ones(size(C,2),1);  % initial condition: input value
options.decim   = decim;
options.ObsEval = ['[posterior0,options0,dim0] = VBA_FinalCheck(posterior,options,dim,suffStat);',...
    'options0.inG.getPhi = 1;',...
    '[gx0] = VBA_evalFun(''g'',posterior0.muX0,posterior0.muPhi,options0.u0,options0,dim0);',...
    'clear posterior0 options0 dim0;'];
% This command is executed after each VB update of the observation
% parameters. It is used to update the persistent variables of the
% generalized observation function.

function [dA,dB,dC,dD] = get_dMatdvec(inF)
if ~isempty(inF.indA)
    dA = dMatdvec(inF.A);
else
    dA = [];
end
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
C = sparse(n,ni);
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
