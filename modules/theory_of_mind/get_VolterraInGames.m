function [out] = get_VolterraInGames(u,tau,param)
% estimate 1st-order Volterra Kernels of agents engaging in games
% FORMAT: [out] = get_VolterraKernels(u,tau,param)
% Volterra kernels can be either estimated in a non-parametric way, or in a
% parametric way (where kernels are assumed to be exponentially decaying
% with lag).
% NB: this function was specifically designed for estimating Volterra
% kernels of agents playing dyadic games with each other. More precisely,
% it estimates the Volterra kernels of the first agent, given both agents'
% actions sequences.
% IN:
%   - u: 2xnt sequence of agents' choices. NB: The Volterra decomposition
%   is performed on u(1,:)!
%   - tau: maximum lag considered
%   - param: flag for parametric (1=default) or non-parametric (0) Volterra
%   decomposition
% OUT: a structure with the following fields:
%   .W: the estimated Volterra kernels
%   .ohm and .A: decay rate and magnitude of Volterra kernels [only for
%   parametric analysis].

try,param;catch,param=1;end

[nu,nt] = size(u);
if param
    g_fname = @g_convSig_param;
else
    g_fname = @g_convSig;
end
inG.dim.n_t = tau;
inG.dim.nu = nu;
dim.n = 0;
dim.n_theta = 0;
dim.n_phi = tau*2 +1;
opt.inG = inG;
opt.priors.muPhi = zeros(dim.n_phi,1);
opt.priors.SigmaPhi = 1e0*eye(dim.n_phi);
opt.sources.type = 1;
opt.DisplayWin = 0;
opt.verbose = 0;

u1 = u(1,1:end-1);
u2 = u(2,1:end-1);
y = u(1,2:end)';
uu = 2*[u1;u2]-1;
%         [uu] = VBA_orth(uu',1)';
[opt.inG.dgdp] = VBA_conv2glm(uu,tau);
[post0,out0] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,opt);

if param
    out.W = extractKernels_param(post0,out0);
    out.ohm = post0.muPhi(1:nu);
    out.A = post0.muPhi(nu+1:2*nu);
else
    out.W = extractKernels(post0,out0);
end




function [sg] = g_convSig_param(x,P,u,in)
% this function evaluates a parametric logistic convolution of the inputs
% function [gx,dgdx,dgdp] = g_convSig_param(x,P,u,in)
% IN:
%   - x: [useless]
%   - P: (K*2+1)x1 vector of kernel parameters.
%   - u: (K*nt)x1 vectorized input to the system
%   - in: [useless]
% OUT:
%   - sg: the predicted system's output (in terms of p(y(t)=1))
% SEE ALSO: g_convSig

tau = in.dim.n_t;
nu = in.dim.nu;
W = zeros(tau*nu,1);
for i=1:nu
    W((i-1)*tau+1:i*tau) = P(i+nu)*exp(-exp(P(i)).*[1:tau]);
end
[sg,dsdx,dsdp] = g_convSig(x,[VBA_vec(W);P(2*nu+1)],u,in);


function [W] = extractKernels_param(posterior,out)
inG = out.options.inG;
tau = inG.dim.n_t;
nu = inG.dim.nu;
W = zeros(tau,nu);
for i=1:nu
    W(:,i) = posterior.muPhi(i+nu)*exp(-exp(posterior.muPhi(i)).*[1:tau]);
end

