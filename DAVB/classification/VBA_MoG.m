function [posterior,out] = VBA_MoG(y,K,options)

% function [posterior,out] = VBA_MoG(y,K,options)
% This function is a VB scheme for the MoG classifier.
% IN:
%   - y : the pxn data matrix, where n is the number of
%   'observations' to be classified, and p is the dimension of each of
%   these observations
%   - K : the maximum number of classes
%   - options : a structure containing infos regarding initialization and
%   inversion scheme (can be left empty -> defaults):
%       .priors: structure containing the following fields:
%           .muEta: pxK prior mean of the components 1st-order moment
%           .SigmaEta: Kx1 cell array of prior pxp covariance matrices of
%           the components 1st-order moment
%           .a_gamma: Kx1 shape parameter of the prior on the components'
%           2nd-order moment
%           .b_gamma: Kx1 scale parameter of the prior on the components'
%           2nd-order moment
%           .d: Kx1 vector of prior counts for each component
%       .TolFun = min increment of the free energy {1e-3}
%       .MaxIter = maximum number of iterations {200}
%       .MinIter = minimum number of iterations {50}
%       .verbose: flag for plotting the clustering results onto the
%       eigenspace of the data
%       .ARD : binary entry that enforces or not the component's death
%       process ({1})
%
% OUT:
%   - posterior: structure containing the following fields:
%       .muEta: pxK prior mean of the components 1st-order moment
%       .SigmaEta: Kx1 cell array of posterior pxp covariance matrices of
%       the components 1st-order moment
%       .a_gamma: Kx1 shape parameter of the posterior on the components'
%       2nd-order moment
%       .a_gamma: Kx1 scale parameter of the posterior on the components'
%       2nd-order moment
%       .d: Kx1 vector of posterior counts for each component
%       .Z: nxK matrix containing the estimated probability, for each
%       voxel, to belong to each of the clusters
%   - out:




% Fill in default options
options.tStart = tic;
[dim.p,dim.n] = size(y);
dim.K = K;
options.dim = dim;
if ~isfield(options,'TolFun')
    options.TolFun = 1e-2;
end
if ~isfield(options,'MaxIter')
    options.MaxIter = 200;
end
if ~isfield(options,'MinIter')
    options.MinIter = 50;
end
if ~isfield(options,'verbose')
    options.verbose = 0;
end
if ~isfield(options,'minSZ')
    options.minSZ = 0e-2/K;
end

% Fill default priors
if ~isfield(options,'priors')
    options.priors = [];
end
priors = options.priors;
if ~isfield(priors,'muEta')
    priors.muEta = zeros(dim.p,dim.K);
end
if ~isfield(priors,'SigmaEta')
    priors.SigmaEta = cell(dim.K,1);
    for k=1:K
        priors.SigmaEta{k} = 1e0*eye(dim.p);
    end
end
if ~isfield(priors,'a_gamma')
    priors.a_gamma = ones(dim.K,1);
end
if ~isfield(priors,'b_gamma')
    priors.b_gamma = ones(dim.K,1);
end
if ~isfield(priors,'d')
    priors.d = ones(dim.K,1);
end
options.priors = priors;

% Initialization
my = mean(y,2);
y = y - repmat(my,1,dim.n);
posterior = priors;
% perturb MoG modes
[eta,Z] = dummyHierarchical(y,K,options);
posterior.muEta = eta;
suffStat.iS0 = cell(dim.K,1);
for k=1:K
    suffStat.iS0{k} = pinv(priors.SigmaEta{k});
end
suffStat.Ed2 = Ed2(y,posterior);

% derive prior on class labels
suffStat.t = zeros(dim.K,dim.n);
posterior.z = zeros(dim.K,dim.n);
Elgam = ElogV(priors);
Egam = EV(priors);
Elr = ElogR(posterior);
suffStat.t = -dim.p*log(2*pi)/2 + repmat(dim.p*Elgam/2+Elr,1,dim.n) - diag(Egam)*suffStat.Ed2/2;
tmp = exp(suffStat.t - repmat(max(suffStat.t,[],1),dim.K,1));
posterior.z = tmp./repmat(sum(tmp,1),dim.K,1);

% for i=1:dim.n
%     suffStat.t(:,i) = -dim.p*log(2*pi)/2 + dim.p*Elgam/2 - Egam.*suffStat.Ed2(:,i)/2 + Elr;
%     tmp = suffStat.t(:,i) - max(suffStat.t(:,i));
%     posterior.z(:,i) = exp(tmp)./sum(exp(tmp));
% end
% derive Free Energy at priors
F = [];%FE(posterior,priors,suffStat,y);

handles.hf = figure('name','Convergence monitoring');
for i=1:9
    handles.ha(i) = subplot(3,3,i,'parent',handles.hf);
end
imagesc(dist(y),'parent',handles.ha(6))
title(handles.ha(6),'data distance')
[u,s,v,] = svd(y,0);
yp = s(1:2,:)*v';
getSubplots

% Main VB scheme
stop = 0;
it = 1;
while ~stop
    
    % display results
    plot(handles.ha(1),F,'ro')
    title(handles.ha(1),'convergence monitoring')
    imagesc(posterior.z,'parent',handles.ha(2))
    title(handles.ha(2),'class labels')
    imagesc(suffStat.Ed2,'parent',handles.ha(3))
    title(handles.ha(3),'distance data - components'' modes')
    imagesc(posterior.muEta,'parent',handles.ha(4))
    title(handles.ha(4),'components'' modes')
    bar(posterior.a_gamma./posterior.b_gamma,'parent',handles.ha(5))
    title(handles.ha(5),'components'' precisions')
    set(handles.ha(7),'nextplot','replace')
    plot(yp(1,:),yp(2,:),'.','parent',handles.ha(7))
    set(handles.ha(7),'nextplot','add')
    mup = u(:,1:2)'*posterior.muEta;
    plot(mup(1,:),mup(2,:),'k+','parent',handles.ha(7))
    drawnow
    pause(0.2)
    % counts for each class
    sumZ = posterior.z*ones(dim.n,1);
    
    % ARD: components' death
    k0 = find(sumZ < options.minSZ);
    if ~isempty(k0) && it > 1
        posterior.a_gamma(k0) = [];
        posterior.b_gamma(k0) = [];
        priors.a_gamma(k0) = [];
        priors.b_gamma(k0) = [];
        priors.muEta(:,k0) = [];
        priors.SigmaEta(k0) = [];
        posterior.muEta(:,k0) = [];
        posterior.SigmaEta(k0) = [];
        posterior.d(k0) = [];
        priors.d(k0) = [];
        posterior.z(k0,:) = [];
        for i=1:dim.n
            posterior.z(:,i) = posterior.z(:,i)./sum(posterior.z(:,i));
        end
        sumZ = posterior.z*ones(dim.n,1);
        suffStat.Ed2(k0,:) = [];
        suffStat.iS0(k0) = [];
        suffStat.t(k0,:) = [];
        dim.K = length(sumZ);
        VBA_disp([num2str(length(k0)),' component''s death: K = ',num2str(dim.K),'.'],options)
    end
    
    % MoG precision parameters
    posterior.a_gamma = priors.a_gamma + dim.p*sumZ/2;
    posterior.b_gamma = priors.b_gamma + sum(posterior.z.*suffStat.Ed2,2);
    
    % MoG class frequencies
    posterior.d = priors.d + sumZ;
    
    % MoG modes
    Egam = EV(posterior);
    for k=1:dim.K
        iS = suffStat.iS0{k} + Egam(k)*sumZ(k)*eye(dim.p);
        posterior.SigmaEta{k} = VB_inv(iS);
        posterior.muEta(:,k) = posterior.SigmaEta{k}*(suffStat.iS0{k}*priors.muEta(:,k)+Egam(k)*y*posterior.z(k,:)');
    end
    suffStat.Ed2 = Ed2(y,posterior);
    
    % MoG class labels
    Elgam = ElogV(posterior);
    Egam = EV(posterior);
    Elr = ElogR(posterior);
    suffStat.t = -dim.p*log(2*pi)/2 + repmat(dim.p*Elgam/2+Elr,1,dim.n) - diag(Egam)*suffStat.Ed2/2;
    tmp = exp(suffStat.t - repmat(max(suffStat.t,[],1),dim.K,1));
    posterior.z = tmp./repmat(sum(tmp,1),dim.K,1);
    
    % Free Energy
    F = [F,FE(posterior,priors,suffStat,y)];
    
    %- Convergence monitoring
    stop = checkStop(it,F,options);
    it = it +1;
    
end

posterior.muEta = posterior.muEta + repmat(my,1,dim.K);

% wraps up the ou structure for display purposes
out.dt = toc(options.tStart);
out.options = options;
out.suffStat = suffStat;
out.date = clock;
out.F = F;
out.it = it;

% SUBFUNCTIONS

function stop = checkStop(it,F,options)
% checks stopping criteria
stop = 0;
if it<options.MinIter
    return
end
dF = F(it) - F(it-1);
if isweird(F) || abs(dF)<=options.TolFun || it>=options.MaxIter
    stop = 1;
end


function F = FE(posterior,priors,suffStat,y)
[p,n] = size(y);
K = size(priors.muEta,2);
Elgam = ElogV(posterior);
Egam = EV(posterior);
Elogr = ElogR(posterior);
F = - n*p*log(2*pi)/2 - p*K*log(2*pi)/2 + p*K*log(2*pi*exp(1))/2 ...
    + gammaln(sum(priors.d)) - sum(gammaln(priors.d)) + priors.d'*Elogr ...
    + priors.a_gamma'*log(priors.b_gamma) - sum(gammaln(priors.a_gamma)) ...
    + (priors.a_gamma-1)'*Elgam - priors.b_gamma'*Egam ...
    + sum(posterior.a_gamma - log(posterior.b_gamma) + gammaln(posterior.a_gamma) + (1-posterior.a_gamma).*psi(posterior.a_gamma)) ...
    - gammaln(sum(posterior.d)) + sum(gammaln(posterior.d)) ...
    - (posterior.d-1)'*Elogr;
for k=1:K
    iS0S = suffStat.iS0{k}*posterior.SigmaEta{k};
    dmu = priors.muEta(:,k) - posterior.muEta(:,k);
    F = F + VBA_logDet(iS0S)/2 - dmu'*suffStat.iS0{k}*dmu/2 + trace(iS0S)/2;
end
F = F + sum(sum(posterior.z.*suffStat.t)) + sum(sum(posterior.z.*log(posterior.z+eps)));


function Elv = ElogV(p)
Elv = psi(p.a_gamma) - log(p.b_gamma);

function Ev = EV(p)
Ev = p.a_gamma./p.b_gamma;

function Elr = ElogR(p)
Elr = psi(p.d) - psi(sum(p.d));


function d2 = Ed2(y,p)
n = size(y,2);
K = size(p.muEta,2);
d2 = zeros(K,n);
trS = zeros(K,1);
for k=1:K
    trS(k) = trace(p.SigmaEta{k});
end
for i=1:n
    tmp = repmat(y(:,i),1,K) - p.muEta;
    for k=1:K
        d2(k,i) = tmp(:,k)'*tmp(:,k) + trS(k);
    end
end

function [flag] = isweird(X)
flag = 0;
if any(isinf(X(:))) || any(isnan(X(:))) || any(~isreal(X(:)))
    flag = 1;
end

function [iy,ix,Max] = maxMat(A)
[my,iyx] = max(A,[],1);
[Max,ix] = max(my);
iy = iyx(ix);

function [eta,Z] = dummyHierarchical(y,K,options)
% This function operates a (dummy) hierarchical clustering algorthm, which
% stops whenever the desired number of clusters is attained. This is used
% for initializing the components' modes.
n = size(y,2);
Z = eye(n);
eta = y;
if options.verbose
    fprintf(1,'Initialing components'' modes (hierarchical clustering)...');
    fprintf(1,'%6.2f %%',0)
end
d = -dist(eta);
for i = 1:n-K
    [i1,i2] = maxMat(d-diag(Inf*ones(size(eta,2),1)));
    newZ = eye(size(eta,2));
    newZ(min([i1,i2]),max([i1,i2])) = 1;
    newZ(max([i1,i2]),:) = [];
    Z = newZ*Z;
    eta = y*Z'*diag(sum(Z,2).^-1);
    tmp = repmat(eta(:,min([i1,i2])),1,size(eta,2)) - eta;
    tmp = -sum(tmp.^2,1);
    d(:,max([i1,i2])) = [];
    d(max([i1,i2]),:) = [];
    d(:,min([i1,i2])) = tmp';
    d(min([i1,i2]),:) = tmp;
    if mod(i,2e1) < 1 && options.verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',100*i/(n-K))
    end
end
if options.verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
end




