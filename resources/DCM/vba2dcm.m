function [DCM] = vba2dcm(posterior,out,DCM,TR)
% maps VBA inversion output to DCM structure
% function [DCM] = vba2dcm(posterior,out,DCM,TR)
% IN:
%   - posterior/out: I/O structures of VBA inversion.
%   - DCM: original DCM structure (see spm_dcm_specify.m). If empty, it is
%   created with default arguments.
%   - TR: this is used only if DCM is left empty...
% OUT:
%   - DCM: DCM structure filled in with the results of VBA inversion

nu = size(out.options.inF.C,2);
nreg = size(out.options.inF.C,1);
y = out.y;
if ~exist('DCM','var') || isempty(DCM) % fill in default info
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
    if out.options.microU
        DCM.U.dt = DCM.Y.dt./out.options.decim;
    else
        DCM.U.dt = DCM.Y.dt;
    end
    DCM.U.u = out.u(out.options.inF.confounds.indu,:)';
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
    DCM.delays = zeros(nreg,1);
end
DCM.M.delays = DCM.delays;
DCM.M.TE = DCM.TE;
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
% include priors
[DCM.M.pE,Pp,DCM.M.pC,Vp] = getABCD(out.options.priors,out.options,0);
try
    DCM.M.pE.transit = out.options.priors.muTheta(out.options.inF.ind2);
    DCM.M.pE.decay = out.options.priors.muTheta(out.options.inF.ind4);
end
DCM.M.x = reshape(out.options.priors.muX0,DCM.M.l,[]);

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
[DCM.Ep,DCM.Pp,DCM.Cp,DCM.Vp] = getABCD(posterior,out.options,0);
try
    DCM.Ep.transit = posterior.muTheta(out.options.inF.ind2);
    DCM.Ep.decay = posterior.muTheta(out.options.inF.ind4);
    % DCM.Ep.epsilon = 0;
end

try
    kernels = out.diagnostics.kernels;
catch
    [kernels] = VBA_getVolterraKernels(posterior,out,ceil(32/TR));
end
DCM.H1 = kernels.g.m;
DCM.K1 = kernels.x.m(out.options.inF.n5,:,:);
% [DCM.H1,DCM.K1] = VBA_getKernels(posterior,out,1);
DCM.M.N = size(DCM.H1,2);
DCM.M.dt = TR;


% fill in inputs if augmented DCM
if isfield(out.options.inF,'augment') ...
        && ~isequal(out.options.inF.augment,0) ...
        && ~isequal(DCM.U.u,out.u')
    nu0 = size(DCM.U.u,1);
    DCM.U.dt = out.options.inF.deltat;
    DCM.U.u = out.u';
    for i=1:out.options.inF.augment.nu
        DCM.U.name{nu0+i} = [out.options.inF.augment.btype,' #',num2str(i)];
        DCM.b(:,:,nu0+i) = zeros(nreg,nreg);
    end
    DCM.c = out.options.inF.C;
end

% finally store VBA output in DCM
DCM.VBA.out = out;
DCM.VBA.posterior = posterior;



function [Ep,Pp,Cp,Vp] = getABCD(posterior,options,t)

inF = options.inF;
Theta = posterior.muTheta;
V = diag(posterior.SigmaTheta);
fullV = posterior.SigmaTheta;
nreg = size(inF.A,1);
nu = size(inF.C,2);
indself = inF.indself;
iself = find(eye(nreg)~=0);
try
%     indhemo = [inF.ind1(:);inF.ind2(:);inF.ind3(:);inF.ind4(:);inF.ind5(:)];
    indhemo = [inF.ind2(:);inF.ind4(:)];
    iHemo = nreg^2*(1+nu+nreg)+nreg*nu+1:nreg^2*(1+nu+nreg)+nreg*(nu+2);
    Cp = zeros(nreg^2*(1+nu+nreg)+nreg*(nu+2)); % 5 hemo params per region!
catch
    indhemo = [];
    iHemo = [];
    Cp = zeros(nreg^2*(1+nu+nreg)+nreg*nu); % 5 hemo params per region!
end

esc = exp(Theta(indself));
vsc = esc^2*V(indself);

Cp(iself,iself) = vsc;
Cp(iHemo,iHemo) = fullV(indhemo,indhemo);
Cp(iself,iHemo) = esc.*repmat(fullV(indself,indhemo),nreg,1);
Cp(iHemo,iself) = esc.*repmat(fullV(indhemo,indself),1,nreg);

Ep.A = zeros(nreg,nreg);
Pp.A = NaN.*ones(nreg,nreg);
Vp.A = zeros(nreg,nreg);
indA = inF.indA;
if ~isempty(indA)
    [ps] = getPPMS(Theta(indA),V(indA),t);
    nzA = find(inF.A~=0);
    Ep.A(nzA) = Theta(indA);
    try
        Pp.A(nzA) = ps;
    end
    Vp.A(nzA) = V(indA);
    Cp(nzA,nzA) = fullV(indA,indA);
    Cp(iself,nzA) = esc.*repmat(fullV(indself,indA),nreg,1);
    Cp(nzA,iself) = esc.*repmat(fullV(indA,indself),1,nreg);
    Cp(iHemo,nzA) = fullV(indhemo,indA);
    Cp(nzA,iHemo) = fullV(indA,indhemo);
end
scI = esc.*eye(nreg);
Ep.A = Ep.A - scI;
if vsc ~= 0
    [psc] = getPPMS(esc,vsc,t);
    Pp.A(scI~=0) = psc;
    Vp.A(scI~=0) = vsc;
else
    Pp.A(scI~=0) = NaN;
end

indB = inF.indB;
offsetB = 0;
for i=1:nu
    Ep.B(:,:,i) = zeros(nreg,nreg);
    Pp.B(:,:,i) = NaN.*ones(nreg,nreg);
    Vp.B(:,:,i) = zeros(nreg,nreg);
    offsetB = offsetB + nreg^2;
    if ~isempty(indB{i})
        Bi = inF.B{i};
        nzB = find(Bi~=0);
        Bi(nzB) = Theta(indB{i});
        Ep.B(:,:,i) = Bi;
        [ps] = getPPMS(Theta(indB{i}),V(indB{i}),t);
        Bi(nzB) = ps;
        Pp.B(:,:,i) = Bi;
        Bi(nzB) = V(indB{i});
        Vp.B(:,:,i) = Bi;
        Cp(nzB+offsetB,nzB+offsetB) = fullV(indB{i},indB{i});
        Cp(iself,nzB+offsetB) = esc.*repmat(fullV(indself,indB{i}),nreg,1);
        Cp(nzB+offsetB,iself) = esc.*repmat(fullV(indB{i},indself),1,nreg);
        Cp(iHemo,nzB+offsetB) = fullV(indhemo,indB{i});
        Cp(nzB+offsetB,iHemo) = fullV(indB{i},indhemo);
        if ~isempty(indA)
            Cp(nzA,nzB+offsetB) = fullV(indA,indB{i});
            Cp(nzB+offsetB,nzA) = fullV(indB{i},indA);
        end
    end
end

Ep.C = zeros(nreg,nu);
Pp.C = NaN.*ones(nreg,nu);
Vp.C = zeros(nreg,nu);
indC = inF.indC;
offsetC = nreg^2*(nu+1);
if ~isempty(indC)
    [ps] = getPPMS(Theta(indC),V(indC),t);
    nzC = find(inF.C~=0);
    Ep.C(nzC) = Theta(indC);   
    Pp.C(nzC) = ps;
    Vp.C(nzC) = V(indC);  
    Cp(nzC+offsetC,nzC+offsetC) = fullV(indC,indC);
    Cp(iself,nzC+offsetC) = esc.*repmat(fullV(indself,indC),nreg,1);
    Cp(nzC+offsetC,iself) = esc.*repmat(fullV(indC,indself),1,nreg);
    Cp(iHemo,nzC+offsetC) = fullV(indhemo,indC);
    Cp(nzC+offsetC,iHemo) = fullV(indC,indhemo);
    if ~isempty(indA)
        Cp(nzA,nzC+offsetC) = fullV(indA,indC);
        Cp(nzC+offsetC,nzA) = fullV(indC,indA);
    end
    offsetB = 0;
    for i=1:nu
        offsetB = offsetB + nreg^2;
        if ~isempty(indB{i}) && ~isempty(find(inF.B{i}~=0))
            nzB = find(inF.B{i}~=0);
            Cp(nzC+offsetC,nzB+offsetB) = fullV(indC,indB{i});
            Cp(nzB+offsetB,nzC+offsetC) = fullV(indB{i},indC);
        end
    end
end

indD = inF.indD;
offsetD = nreg^2*nu + nreg*nu;
for i=1:nreg
    Ep.D(:,:,i) = zeros(nreg,nreg);
    Pp.D(:,:,i) = NaN.*ones(nreg,nreg);
    Vp.D(:,:,i) = zeros(nreg,nreg);
    offsetD = offsetD + nreg^2;
    if ~isempty(indD{i})
        Di = inF.D{i};
        nzD = find(Di~=0);
        Di(nzD) = Theta(indD{i});
        Ep.D(:,:,i) = Di;
        [ps] = getPPMS(Theta(indD{i}),V(indD{i}),t);
        Di(nzD) = ps;
        Pp.D(:,:,i) = Di;
        Di(nzD) = V(indD{i});
        Vp.D(:,:,i) = Di;
        Cp(nzD+offsetD,nzD+offsetD) = fullV(indD{i},indD{i});
        Cp(iself,nzD+offsetD) = esc.*repmat(fullV(indself,indD{i}),nreg,1);
        Cp(nzD+offsetD,iself) = esc.*repmat(fullV(indD{i},indself),1,nreg);
        Cp(iHemo,nzD+offsetD) = fullV(indhemo,indD{i});
        Cp(nzD+offsetD,iHemo) = fullV(indD{i},indhemo);
        if ~isempty(indA)
            Cp(nzA,nzD+offsetD) = fullV(indA,indD{i});
            Cp(nzD+offsetD,nzA) = fullV(indD{i},indA);
        end
        offsetB = 0;
        for j=1:nu
            offsetB = offsetB + nreg^2;
            if ~isempty(indB{j}) && ~isempty(find(inF.B{j}~=0))
                nzB = find(inF.B{j}~=0);
                Cp(nzD+offsetD,nzB+offsetB) = fullV(indD{i},indB{j});
                Cp(nzB+offsetB,nzD+offsetD) = fullV(indB{j},indD{i});
            end
        end
        offsetC = nreg^2*(nu+1);
        if ~isempty(indC) && ~isempty(find(inF.C~=0))
            nzC = find(inF.C~=0);
            Cp(nzC+offsetC,nzD+offsetD) = fullV(indC,indD{i});
            Cp(nzD+offsetD,nzC+offsetC) = fullV(indD{i},indC);
        end
    end
end


function [ps] = getPPMS(m,v,t)
n = length(m);
ps = zeros(n,1);
for i=1:n
    [ps(i)] = VBA_PPM(abs(m(i)),v(i),t);
end
