function [m,s,w] = getMoG4N01(nmog,cmog,verbose)
% fits a MoG model to the N(0,1) density
% function [m,s,w] = getMoG4N01(nmog,verbose)
% Let us posit the following mixture of gaussians (MoG) model:
%   N(0,1) = sum_n w(n) N(m(n),s(n)) + e
% where the sum is for all n smaller than nmog, which is the number of
% gaussian components in the MoG model.
% This function inverts the above MoG model under he following constraints:
% - the variances s(n) are the same for all components
% - the components are paired, i.e. for any component with mean m(n), there
% is one with mean -m(n), and both have the same weight. NB: when nmog is
% odd, there is a central component with mean fixed to zero.
% - each (paired) mean is forced to be strictly greater than the preceeding
% one, i.e.: m(n+1) > m(n).
% IN:
%   - nmog: number of univariate gaussian components in the MoG model
%   - cmog: the class of constained MoG model
%   - verbose: falg for verbose mode
% OUT:
%   - m: nmogx1 vector of means
%   - s: nmogx1 vector of variances
%   - w: nmogx1 vector of weights

try, cmog; catch, cmog = 1; end
try, verbose; catch, verbose = 0; end

% trivial input
if isequal(nmog,1)
    m = 0;
    s = 1;
    w = 1;
    return
end

% evaluate N(0,1) on the grid
inG.gri = -20:2e-2:20;  % grid on which the gbf is evaluated
gx = exp(-0.5*(inG.gri).^2)./sqrt(2*pi);
y = gx(:);% + sqrt(1e-5).*randn(length(gx),1);

% set up VB inversion routine
g_fname = @constrainedMoG; % observation function
if isequal(nmog/2,floor(nmog/2))
    % even #components
    np = nmog/2;
    if cmog == 1
        n_theta = 2*np+1;
    elseif cmog == 2
        n_theta = 2+np;
    end
else
    % odd #components
    np = (nmog-1)./2;
    if cmog == 1
        n_theta = 2*np+2;
        inG.q0 = 2+2*np;
    elseif cmog == 2
        n_theta = 2+np;
    end
end
inG.s = 1;
inG.mu = 2:1+np;
if cmog == 1
    inG.q = 2+np:1+2*np;
elseif cmog == 2
    inG.q = 2+np;
end
inG.nmog = nmog;
inG.cmog = cmog;
priors.muPhi = zeros(n_theta,1);
priors.SigmaPhi = 1e-2*speye(n_theta);
if nmog==2 || cmog == 2
    % if nmog=2, then weights=0.5 ...
    priors.SigmaPhi(inG.q,inG.q) = 0;
end
priors.a_sigma = 1e5;
priors.b_sigma = 1;
options.priors = priors;
options.inG = inG;
options.GnFigs = 0;
options.verbose = 0;
options.DisplayWin = verbose;
dim.n_phi = n_theta;
dim.n_theta = 0;
dim.n=0;

% Call inversion routine
[p,o] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);


% get MoG sufficient statistics
Phi = p.muPhi;
mu = cumsum(exp(Phi(inG.mu)));
s = sigm(Phi(inG.s))*ones(nmog,1);
w = zeros(nmog,1);
m = zeros(nmog,1);
if isequal(nmog/2,floor(nmog/2))
    % even #components
    if cmog ==1
        q = exp(Phi(inG.q));
    elseif cmog == 2
        q = exp(-0.5*mu.^2.*exp(Phi(inG.q)));
    end
    q = q./2.*sum(q);
    w(1:2:nmog-1) = q;
    w(2:2:nmog) = q;
    m(1:2:nmog-1) = mu;
    m(2:2:nmog) = -mu;
else
    % odd #components
    if cmog == 1
        q = exp(Phi(inG.q));
        q0 = exp(Phi(inG.q0));
    elseif cmog == 2
        q = exp(-0.5*mu.^2.*exp(Phi(inG.q)));
        q0 = 1;
    end
    q = q./(2*sum(q)+q0);
    q0 = q0./(2*sum(q)+q0);
    w(1) = q0;
    w(2:2:nmog-1) = q;
    w(3:2:nmog) = q;
    m(1) = 0;
    m(2:2:nmog-1) = mu;
    m(3:2:nmog) = -mu;
end

% display fit?
if verbose
    ha = plot_MoG(m,s,w,inG);
    plot(ha,inG.gri,gx,'g--')
    getSubplots
end

