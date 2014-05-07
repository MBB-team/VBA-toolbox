function [H1,K1,tgrid] = VBA_getKernels(posterior,out,dcm)
% Dummy derivation of the system's response kernels
% function [H1,K1,tgrid] = VBA_getKernels(posterior,out,dcm)
% IN:
%   - posterior,out: the output of the system inversion
%   - dcm: flag for dcm (does not compute hemodynamic states kernels)
% OUT:
%   - H1: the (pxtxnu) output impulse response function, where p is the
%   dimension of the data, t is the number of time samples and nu is the
%   number of inputs to the system
%   - K1: the (nxtxnu) state impulse response function, where n is the
%   number of states. NB: for DCM models, this is the neural impulse
%   response function...
%   - tgrid: the time grid over which the kernels are estimated

if isequal(out.dim.n_t,1) || out.dim.u < 1 || isempty(out.options.f_fname)
    % not a dynamical system
    H1 = [];
    K1 = [];
    tgrid = [];
    return
end

if nargin <3
    dcm = 0;
else
    dcm = ~~dcm;
end

nu = out.dim.u;

n = out.dim.n;
p = out.dim.p;
out.options.microU = 1; % for impulse response functions...
if dcm
    % remove confounds,...
    out.options.inG.confounds.X0 = [];
    % only look at neuronal states,...
    n = p; %size(out.options.inF.C,1);
    % and get kernels over 16 secs
    TR = out.options.inF.deltat*out.options.decim;
    out.options.dim.n_t = ceil(16./TR);
end
nt = out.options.dim.n_t*out.options.decim;

% ensure steady state initial conditions and throw away state noise
% estimate
posterior.muX0 = zeros(size(posterior.muX0));
out.suffStat.dx = [];

% pre-allocate response kernels
H1 = zeros(p,nt,nu);
K1 = zeros(n,nt,nu);

% derive kernels by integrating the system
gotit = 0;
for i=1:nu
    try
        U = zeros(nu,nt);
        U(i,1) = 1;
        % get output impulse response
        [x,gx,tgrid] = VBA_microTime(posterior,U,out);
        H1(:,:,i) = gx(:,2:end);
        if dcm && isfield(out.options.inF,'n5')
            K1(:,:,i) = x(out.options.inF.n5,2:end);
        else
            K1(:,:,i) = x(:,2:end);
        end
        gotit = 1;
    end
end
% clean up kernels
if ~gotit
    H1 = [];
    K1 = [];
    tgrid = [];
else
    tgrid = tgrid(2:end);
    K1(abs(K1)<=1e-8) = 0;
    H1(abs(H1)<=1e-8) = 0;
end

