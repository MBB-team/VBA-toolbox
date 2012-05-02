function [muy,Vy,m,V,w] = splitLaplace(u,f_fname,g_fname,dim,options,nmog,nosplit)
% returns the split-Laplace approximation to the prior predictive density
% function [muy,Vy] = splitLaplace(u,f_fname,g_fname,dim,options,nmog)
% IN:
%   - u: experimentally controlled input (design)
%   - f_fname: the evolution function
%   - g_fname: the observation function
%   - dim: the model dimension structure
%   - options: the options structure
%   - nmog: number of MoG components in the split-Laplace approx {2}
%   - nosplit: structure containing the following fields:
%       .phi: vector of indices of the observation prameters for which no
%       split is performed
%       .theta: vector of indices of the evolution prameters for which no
%       split is performed
%       .x0: vector of indices of the initial conditions for which no split
%       is performed
% OUT:
%   - muy: the 1st-order moment of the prior predictive density.
%   - Vy: the second-order moment of the prior predictive density.
% SEE ALSO: getLaplace


% 0- check input and dimensions to split
try; nmog; catch; nmog = 2; end
options.priors.a_alpha = 0; % to bypass ODE transform in VBA_check.m
options.verbose = 0; % to quicken VBA_check.m
[options,u,dim] = VBA_check([],u,f_fname,g_fname,dim,options);
u = VBA_getU(u,options,dim,'back2micro');

try
    nosplit.phi = unique([nosplit.phi(:);...
        find(diag(options.priors.SigmaPhi)==0)]);
catch
    nosplit.phi = find(diag(options.priors.SigmaPhi)==0);
end
try
    nosplit.theta = unique([nosplit.theta(:);...
        find(diag(options.priors.SigmaTheta)==0)]);
catch
    nosplit.theta = find(diag(options.priors.SigmaTheta)==0);
end
try
    nosplit.x0 = unique([nosplit.x0(:);...
        find(diag(options.priors.SigmaX0)==0)]);
catch
    nosplit.x0 = find(diag(options.priors.SigmaX0)==0);
end
out = [...
    nosplit.phi
    nosplit.theta+dim.n_phi
    nosplit.x0+dim.n_phi+dim.n_theta
    ];
np = dim.n_phi+dim.n_theta+dim.n;
in = setdiff([1:np]',out);
if isempty(in)
    [muy,Vy] = getLaplace(u,f_fname,g_fname,dim,options);
    m = muy;
    V = {Vy};
    w = 1;
    return
end
    
% 1- get compact prior mean and covariance matrix
Sigma = zeros(dim.n_phi+dim.n_theta+dim.n,dim.n_phi+dim.n_theta+dim.n);
Mu = zeros(dim.n_phi+dim.n_theta+dim.n,1);
if dim.n_phi > 0
    Sigma(1:dim.n_phi,1:dim.n_phi) = options.priors.SigmaPhi;
    Mu(1:dim.n_phi) = options.priors.muPhi;
end
if dim.n_theta > 0
    Sigma(dim.n_phi+1:dim.n_phi+dim.n_theta,...
        dim.n_phi+1:dim.n_phi+dim.n_theta) = options.priors.SigmaTheta;
    Mu(dim.n_phi+1:dim.n_phi+dim.n_theta) = options.priors.muTheta;
end
if dim.n > 0
    Sigma(dim.n_phi+dim.n_theta+1:end,dim.n_phi+dim.n_theta+1:end) = ...
        options.priors.SigmaX0;
    Mu(dim.n_phi+dim.n_theta+1:end) = options.priors.muX0;
end
% get square-root of full covariance matrix
sqrtS = getISqrtMat(Sigma(in,in),0);

% 2- fit MoG (with nmog components) to the N(0,1) density
[m0,s0,w0] = getMoG4N01(nmog,1,0);

% 3- get n-draws from nmog-urn
[C] = get_nkdraws(nmog,length(in),0);

% 4- loop over all possible n-draws
nd = size(C,2);
w = zeros(1,nd);
m = zeros(dim.p*dim.n_t,nd);
V = cell(nd,1);
Mu0 = Mu;
Sigma0 = Sigma;
fprintf(1,'%6.2f  %%',0)
for i=1:nd
    mi = m0(C(:,i));
    Si = diag(s0(C(:,i)));
    w(i) = prod(w0(C(:,i)));
    Mu0(in) = sqrtS*mi + Mu(in);
    Sigma0(in,in) = sqrtS*Si*sqrtS';
    [options.priors] = getPriors4split(Mu0,Sigma0,dim,options.priors);
    [m(:,i),V{i}] = getLaplace(u,f_fname,g_fname,dim,options);
    fprintf(1,'\b\b\b\b\b\b\b\b')
    fprintf(1,'%6.2f %%',100*i/nd)
end

% 5- match moments of the MoG
w = w./sum(w);
muy = sum(m*diag(w(:)),2);
Vy = zeros(dim.p*dim.n_t,dim.p*dim.n_t);
for i=1:size(C,2)
    dm = m(:,i) - muy;
    Vy = Vy + w(i).*(dm*dm'+V{i});
end


function [priors] = getPriors4split(Mu0,Sigma0,dim,priors)
if dim.n_phi > 0
    priors.muPhi = Mu0(1:dim.n_phi);
    priors.SigmaPhi = Sigma0(1:dim.n_phi,1:dim.n_phi);
end
if dim.n_theta > 0
    priors.muTheta = Mu0(dim.n_phi+1:dim.n_phi+dim.n_theta);
    priors.SigmaTheta = Sigma0(...
        dim.n_phi+1:dim.n_phi+dim.n_theta,...
        dim.n_phi+1:dim.n_phi+dim.n_theta);
end
if dim.n > 0
    priors.muX0 = Mu0(dim.n_phi+dim.n_theta+1:end);
    priors.SigmaX0 = Sigma0(...
        dim.n_phi+dim.n_theta+1:end,...
        dim.n_phi+dim.n_theta+1:end);
end


