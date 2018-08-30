function [y,u,f_fname,g_fname,dim,options] = dcm2vba(DCM,stochastic,augment,confounds)
% maps DCM input structure to VBA inversion scheme input
% function [y,u,f_fname,g_fname,dim,options] = dcm2vba(DCM,stochastic,augment)
% IN:
%   - DCM: DCM structure (see spm_dcm_specify.m)
%   - stochastic: binary variables that flags the use of either stochastic
%   (1) or deterministic DCM (0)
%   - augment: a structure that serves to buid input basis function sets,
%   with the following fields:
%       .btype: either 'Fourier' or 'RBF'. This specifies the type of input
%       basis function set
%       .nu: the number of basis functions
% OUT:
%   - y: the fMRI data, adjusted for confounds and delays
%   - u: the system's inputs, in microtime resolution. Note: input basis
%   functions are concatenated with controlled inputs in u.
%   - f_fname: handle of the evolution function for VBA inversion
%   - g_fname: handle of the observation function for VBA inversion
%   - dim: model's dimension structure
%   - options: I/O structure for VBA inversion

try; stochastic; catch, stochastic=0; end
try; augment; catch, augment=0; end
try; confounds; catch, confounds=0; end

% hard coded parameters
microDT = 1e-1; % micro-time resolution (in sec)
reduced_f = 0;%1; % do not fix HRF parameters


% Unpack DCM specification and derive options structure
[ny,nreg] = size(DCM.Y.y);
dt = DCM.U.dt;
TR = DCM.Y.dt;
uu = DCM.U.u';
[nu,nt] = size(uu);
A = DCM.a - eye(nreg);
if nu > 0
    for i=1:nu
        B{i} = DCM.b(:,:,i);
    end
else
    B{1} = zeros(nreg,nreg);
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


if ~isempty(augment) && ~isequal(augment,0)
    % add input basis function sets
    if nu==0
        nt = TR*ny/dt;
        uu = [];
        C = [];
        B = cell(0);
    end
    xu = get_U_basis(nt*dt,dt,augment.nu,augment.btype);
    uu = [uu;xu];
    C = [C,ones(nreg,augment.nu)];
    for i=length(B)+1:size(C,2)
        B{i} = zeros(nreg,nreg);
    end
    nu = size(uu,1);
end


% Build priors and finalize options and dim structures
f_fname = @f_DCMwHRF; % DCM evolution function
g_fname = @g_HRF3; % DCM observation function
[options] = prepare_fullDCM(A,B,C,D,TR,microDT,0);
options.priors = getPriors(nreg,ny,options,reduced_f,stochastic);
options.DisplayWin = 0; % VBA inversion window is shown
options.backwardLag = ceil(16/TR); % spans 16secs of fMRI data
options.init0 = 0; % use dDCM parameter estimate for sDCM inversion
options.microU = 1; % inputs are defined at microtime resolution
options.inF.augment = augment;
options.inF.linearized = 1; % no nonlinearity in the HRF balloon model
dim.n_theta = options.inF.ind5(end);
dim.n_phi = options.inG.ind2(end);
dim.n = 5*nreg;
dim.p = nreg;
dim.n_t = ny;



% Resample inputs on microtime integration grid
[u] = resampleU(uu,dt,nt,nu,ny,options);
u = [zeros(nu,1),u];
if size(u,1)==0
    u = zeros(1,size(u,2));
end


% Get data and add in confounds
y0 = DCM.Y.y';
try
    X0 = DCM.Y.X0; % minimal confounds: slow drifts
catch
    X0 = [];
end
if ~isequal(confounds,0)
    if isequal(confounds,'find')
        [X0,isYout] = VBA_spm_Xadjust(DCM.xY(1).SPMfile,DCM.xY(1).VOIfile);
    else
        X0 = confounds.X0;
        isYout = confounds.isYout;
    end
    for t=1:length(isYout)
        options.priors.iQy{isYout(t)} = 0.*options.priors.iQy{isYout(t)};
    end
    X0 = [X0,DCM.Y.X0];
end
[y0] = adjustY(y0,X0);
[u,options,dim] = addConfounds2dcm(X0,u,options,dim);

% Deal with delays
y = zeros(size(y0));
for i=1:nreg
    try
        D = round(DCM.delays(i)./TR);
        y(i,:) = [zeros(1,D),y0(i,1:end-D)];
    catch % no delay specified
        y(i,:) = y0(i,:);
    end
end




function [u] = resampleU(uu,dt,nt,nu,ny,options)
microDT = options.inF.deltat;
grid1 = 0:dt:dt*(nt-1);
grid2 = 0:microDT:microDT*options.decim*ny;
u = zeros(nu,length(grid2)-1);
for i=1:length(grid2)-1
    [tmp,ind1] = min(abs(grid2(i)-grid1));
    [tmp,ind2] = min(abs(grid2(i+1)-grid1));
    u(:,i) = mean(uu(:,ind1:ind2),2);
end
% [u,alpha] = spm_resample(full(uu),dt/microDT);

function [y] = adjustY(y,X0)
if ~isempty(X0)
    iX0 = pinv(X0'*X0)*X0';
    for i=1:size(y,1)
        beta = iX0*y(i,:)';
        yc = X0*beta;
        y(i,:) = y(i,:) - yc';
    end
end


