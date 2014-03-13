function [DCM] = exportDCMfromVBNLSS(posterior,out,DCM,TR)

% This function fills in a DCM structure with the output of VBNLLS
% (i.e. stochastic DCM) inversion

nu = size(out.options.inF.C,2);
nreg = size(out.options.inF.C,1);
y = out.y;
if isempty(DCM) % fill in default info
    try
        DCM.Y.dt = TR;
    catch
        DCM.Y.dt = 1;
    end
    try
        DCM.TE = out.options.inG.TE;
    catch
        DCM.TE = 0.04;
    end
    DCM.Y.y = y';
    DCM.U.dt = 1;
    DCM.U.u = out.u';
    for i=1:nreg
        DCM.Y.name{i} = ['region #',num2str(i)];
    end
    for i=1:nu
        DCM.U.name{i} = ['factor #',num2str(i)];
    end
    DCM.a = out.options.inF.A;
    for i=1:nu
        if ~isempty(out.options.inF.indB{i})
            DCM.b(:,:,i) = out.options.inF.B{i};
        else
            DCM.b(:,:,i) = zeros(nreg,nreg);
        end
    end
    DCM.c = out.options.inF.C;
    for i=1:nreg
        if ~isempty(out.options.inF.indD{i})
            DCM.d(:,:,i) = out.options.inF.D{i};
        else
            DCM.d(:,:,i) = zeros(nreg,nreg);
        end
    end
end
DCM.M.m = nu; % number of inputs
DCM.M.n = out.dim.n; % number of hidden states
DCM.M.l = nreg; % number of regions
DCM.M.IS = '[in built Euler integration scheme]';
if isa(out.options.g_fname,'function_handle')
    DCM.M.g = func2str(out.options.g_fname);
else
    DCM.M.g = out.options.g_fname;
end
if isa(out.options.f_fname,'function_handle')
    DCM.M.f = func2str(out.options.f_fname);
else
    DCM.M.f = out.options.f_fname;
end
DCM.options.two_state = 0;
if isinf(out.options.priors.a_alpha) && ...
        isequal(out.options.priors.b_alpha,0)
    sdcm = 0;
else
    sdcm = 1;
end
DCM.options.stochastic = sdcm;
DCM.y = out.suffStat.gx'; % predicted data
DCM.R = y' - DCM.y; % residuals of the model
DCM.F = out.F;

% gets A, B, C and D expectations and PPMs
[A,B0,C,D0,pA,pB0,pC,pD0,VAR] = getABCD(posterior,out,0);
for i=1:nu
    B(:,:,i) = B0{i};
    pB(:,:,i) = pB0{i};
end
for i=1:nreg
    D(:,:,i) = D0{i};
    pD(:,:,i) = pD0{i};
end

DCM.Cp = VAR;

DCM.Ep.A = A;
DCM.Ep.B = B;
DCM.Ep.C = C;
DCM.Ep.D = D;

DCM.Pp.A = pA;
DCM.Pp.B = pB;
DCM.Pp.C = pC;
DCM.Pp.D = pD;

[DCM.H1,DCM.K1] = getKernels(posterior,out,1);
DCM.M.N = size(DCM.H1,2);
DCM.M.dt = out.options.inF.deltat;

% finally store VBA outout structures in DCM
DCM.VBA = struct('out',out,'posterior',posterior);




function [A,B,C,D,pA,pB,pC,pD,VAR] = getABCD(posterior,out,t)

inF = out.options.inF;
Theta = posterior.muTheta;
V = diag(posterior.SigmaTheta);
fullV = posterior.SigmaTheta;
nreg = size(inF.A,1);
nu = size(inF.C,2);
indself = inF.indself;
iself = find(eye(nreg)~=0);
try
    indhemo = [...
        inF.ind1(:)...
        inF.ind2(:)...
        inF.ind3(:)...
        inF.ind4(:)...
        inF.ind5(:)...
        ];
    iHemo = nreg^2*(1+nu+nreg)+nreg*nu+1:nreg^2*(1+nu+nreg)+nreg*(nu+5);
catch
    indhemo = [];
    iHemo = [];
end
VAR = zeros(nreg^2*(1+nu+nreg)+nreg*(nu+5)); % 5 hemo params per region!

esc = exp(Theta(indself));
vsc = esc.^2.*V(indself);

VAR(iself,iself) = vsc;
VAR(iHemo,iHemo) = fullV(indhemo,indhemo);
VAR(iself,iHemo) = esc.*repmat(fullV(indself,indhemo),nreg,1);
VAR(iHemo,iself) = esc.*repmat(fullV(indhemo,indself),1,nreg);

A = inF.A;
pA = NaN.*ones(nreg,nreg);
indA = inF.indA;
if ~isempty(indA)
    [ps] = getPPMS(Theta(indA),V(indA),t);
    nzA = find(A~=0);
    A(nzA) = Theta(indA);
    pA(nzA) = ps;
    VAR(nzA,nzA) = fullV(indA,indA);
    VAR(iself,nzA) = esc.*repmat(fullV(indself,indA),nreg,1);
    VAR(nzA,iself) = esc.*repmat(fullV(indA,indself),1,nreg);
    VAR(iHemo,nzA) = fullV(indhemo,indA);
    VAR(nzA,iHemo) = fullV(indA,indhemo);
end
scI = esc.*eye(nreg);
A = A - scI;
if vsc ~= 0
    [psc] = getPPMS(esc,vsc,t);
    pA(scI~=0) = psc;
else
    pA(scI~=0) = NaN;
end

B = inF.B;
pB = cell(nu,1);
indB = inF.indB;
offsetB = 0;
for i=1:nu
    offsetB = offsetB + nreg^2;
    if ~isempty(indB{i})
        [ps] = getPPMS(Theta(indB{i}),V(indB{i}),t);
        nzB = find(B{i}~=0);
        B{i}(nzB) = Theta(indB{i});
        pB{i} = NaN.*ones(nreg,nreg);
        pB{i}(nzB) = ps;
        VAR(nzB+offsetB,nzB+offsetB) = fullV(indB{i},indB{i});
        VAR(iself,nzB+offsetB) = ...
            esc.*repmat(fullV(indself,indB{i}),nreg,1);
        VAR(nzB+offsetB,iself) = ...
            esc.*repmat(fullV(indB{i},indself),1,nreg);
        VAR(iHemo,nzB+offsetB) = fullV(indhemo,indB{i});
        VAR(nzB+offsetB,iHemo) = fullV(indB{i},indhemo);
        if ~isempty(indA)
            VAR(nzA,nzB+offsetB) = fullV(indA,indB{i});
            VAR(nzB+offsetB,nzA) = fullV(indB{i},indA);
        end
    else
        B{i} = zeros(nreg,nreg);
        pB{i} = NaN.*ones(nreg,nreg);
    end
end

C = inF.C;
pC = NaN.*ones(nreg,nu);
indC = inF.indC;
offsetC = nreg^2*(nu+1);
if ~isempty(indC)
    [ps] = getPPMS(Theta(indC),V(indC),t);
    nzC = find(C~=0);
    C(nzC) = Theta(indC);   
    pC(nzC) = ps;
    VAR(nzC+offsetC,nzC+offsetC) = fullV(indC,indC);
    VAR(iself,nzC+offsetC) = esc.*repmat(fullV(indself,indC),nreg,1);
    VAR(nzC+offsetC,iself) = esc.*repmat(fullV(indC,indself),1,nreg);
    VAR(iHemo,nzC+offsetC) = fullV(indhemo,indC);
    VAR(nzC+offsetC,iHemo) = fullV(indC,indhemo);
    if ~isempty(indA)
        VAR(nzA,nzC+offsetC) = fullV(indA,indC);
        VAR(nzC+offsetC,nzA) = fullV(indC,indA);
    end
    offsetB = 0;
    for i=1:nu
        offsetB = offsetB + nreg^2;
        if ~isempty(indB{i})
            nzB = find(B{i}~=0);
            VAR(nzC+offsetC,nzB+offsetB) = fullV(indC,indB{i});
            VAR(nzB+offsetB,nzC+offsetC) = fullV(indB{i},indC);
        end
    end
end

D = inF.D;
pD = cell(nreg,1);
indD = inF.indD;
offsetD = nreg^2*nu + nreg*nu;
for i=1:nreg
    offsetD = offsetD + nreg^2;
    if ~isempty(indD{i})
        nzD = find(D{i}~=0);
        [ps] = getPPMS(Theta(indD{i}),V(indD{i}),t);
        D{i}(nzD) = Theta(indD{i});
        pD{i} = NaN.*ones(nreg,nreg);
        pD{i}(nzD) = ps;
        VAR(nzD+offsetD,nzD+offsetD) = fullV(indD{i},indD{i});
        VAR(iself,nzD+offsetD) = ...
            esc.*repmat(fullV(indself,indD{i}),nreg,1);
        VAR(nzD+offsetD,iself) = ...
            esc.*repmat(fullV(indD{i},indself),1,nreg);
        VAR(iHemo,nzD+offsetD) = fullV(indhemo,indD{i});
        VAR(nzD+offsetD,iHemo) = fullV(indD{i},indhemo);
        if ~isempty(indA)
            VAR(nzA,nzD+offsetD) = fullV(indA,indD{i});
            VAR(nzD+offsetD,nzA) = fullV(indD{i},indA);
        end
        offsetB = 0;
        for j=1:nu
            offsetB = offsetB + nreg^2;
            if ~isempty(indB{j})
                nzB = find(B{j}~=0);
                VAR(nzD+offsetD,nzB+offsetB) = ...
                    fullV(indD{i},indB{j});
                VAR(nzB+offsetB,nzD+offsetD) = ...
                    fullV(indB{j},indD{i});
            end
        end
        offsetC = nreg^2*(nu+1);
        if ~isempty(indC)
            VAR(nzC+offsetC,nzD+offsetD) = fullV(indC,indD{i});
            VAR(nzD+offsetD,nzC+offsetC) = fullV(indD{i},indC);
        end
    else
        D{i} = zeros(nreg,nreg);
        pD{i} = NaN.*ones(nreg,nreg);
    end
end


function [ps] = getPPMS(m,v,t)
n = length(m);
ps = zeros(n,1);
for i=1:n
    [ps(i)] = VBA_PPM(abs(m(i)),v(i),t,0);
end
