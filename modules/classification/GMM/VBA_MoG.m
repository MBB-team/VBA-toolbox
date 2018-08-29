function [posterior,out] = VBA_MoG(y,K,options)
% function [posterior,out] = VBA_MoG(y,K,options)
% This function is a VB scheme for the MoG classifier.
% IN:
%   - y : the pxn data matrix, where n is the number of 'observations' to be
%   classified, and p is the dimension of each of these observations
%   - K : the maximum number of classes
%   - options : a structure containing infos regarding initialization and
%   inversion scheme (can be left empty -> defaults):
%       .priors: structure containing the following fields:
%           .muEta: pxK prior mean of the components 1st-order moment
%           .SigmaEta: Kx1 cell array of prior pxp covariance matrices of
%           the components 1st-order moment
%           .a_gamma: Kx1 vector of shape parameters for the prior on the
%           components' 2nd-order moment
%           .b_gamma: Kx1 vector of scale parameters for the prior on the
%           components' 2nd-order moment
%           .d: Kx1 vector of prior counts for each component
%       .normalize: flag for data normalization, i.e. remove mean along
%       each dimension, and divide by global standard deviation {1}
%       .TolFun: min increment of the free energy {1e-2}
%       .MaxIter: maximum number of iterations {256}
%       .MinIter: minimum number of iterations {0}
%       .verbose: flag for verbose mode
%       .DisplayWin: flag for plotting the clustering results onto the
%       eigenspace of the data
%       .init: {'hierarchical'}, 'prior' or 'rand'. This is a flag for how
%       to initialize the posterior mean over MoG modes. The dafult options
%       uses a (dummy) hierarchical clustering technique, 'prior' is just
%       the prior mean and 'rand' is a sample under the prior density
%       (useful for multistart variants)
%       .minSumZ: the threshold on the sum of class counts, below which
%       components are removed from the model {1}. To eschew this ARD
%       process, one simply sets options.minSumZ = 0.
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
%       voxel, of belonging to each of the clusters
%   - out:
%       .dt: elapsed time (in sec)
%       .dim: final dimension of the model
%       .options: (filled in) options for MoG VB model inversion
%       .suffStat: internal sufficient statistics
%       .date: date (vector format)
%       .F: history of Free Energy across VB iterations
%       .it: # VB iterations until convergence
%       .normalize: a structure containing the necessary information for
%       de-normalization of the inferred class patterns, i.e.:
%           muEta <-- diag(sqrt(diag(Q)))*muEta + repmat(m,1,dim.K);
%           b_gamma(k) <-- = b_gamma(k).*Q;
%       NB: this is performed automatically by VBA_projectMoG.m, when
%       visualizing the decision boundaries.

% First of all, deal with data normalization
if ~isfield(options,'normalize')
    options.normalize = 0;
end
y0 = y;
if options.normalize
    my = mean(y,2);
    y = y - repmat(my,1,size(y,2));
    sy = std(y,[],2);
    if any(sy==0)
        i0 = find(sy==0);
        y(i0,:) = [];
        y0(i0,:) = [];
        my(i0) = [];
        sy(i0) = [];
        disp('Warning: removing 0 variance data dimensions!')
    else
        i0 = [];
    end
    y = diag(1./sy)*y;
else
    i0 = [];
end

% Fill in default options
options.tStart = tic;
[dim.p,dim.n] = size(y);
dim.K = K;
options.dim = dim;
if ~isfield(options,'TolFun')
    options.TolFun = 1e-2;
end
if ~isfield(options,'MaxIter')
    options.MaxIter = 256;
end
if ~isfield(options,'MinIter')
    options.MinIter = 0;
end
if ~isfield(options,'verbose')
    options.verbose = 1;
end
if ~isfield(options,'DisplayWin')
    options.DisplayWin = 1;
end
if ~isfield(options,'minSumZ')
    options.minSumZ = 1;%1/K;
end
if ~isfield(options,'init')
    options.init = 'hierarchical';
end

% Fill in default priors
if ~isfield(options,'priors')
    options.priors = [];
end
priors = options.priors;
if ~isfield(priors,'muEta')
    priors.muEta = zeros(dim.p,dim.K);
end
priors.muEta(i0,:) = [];
if ~isfield(priors,'SigmaEta')
    priors.SigmaEta = cell(dim.K,1);
    for k=1:K
        priors.SigmaEta{k} = 1e0*eye(dim.p);
    end
end
suffStat.iS0 = cell(dim.K,1);
for k=1:K
    priors.SigmaEta{k}(i0,:) = [];
    priors.SigmaEta{k}(:,i0) = [];
    suffStat.iS0{k} = VBA_inv(priors.SigmaEta{k});
end
if ~isfield(priors,'a_gamma')
    priors.a_gamma = 1e0*ones(dim.K,1);
end
if ~isfield(priors,'b_gamma')
    priors.b_gamma = ones(dim.K,1);
end
if ~isfield(priors,'d')
    priors.d = ones(dim.K,1);
end
options.priors = priors;

if options.verbose
    fprintf('\n')
    fprintf('---')
    fprintf('\n')
end

% Initialization
posterior = priors;
switch options.init
    case 'hierarchical'
        [posterior.muEta,Z] = dummyHierarchical(y,K,options);
    case 'prior'
        posterior.muEta = priors.muEta + eps*randn(dim.p,dim.K);
    case 'rand'
        for k=1:K
            S = VBA_sqrtm (priors.SigmaEta{k});
            posterior.muEta(:,k) = priors.muEta(:,k) + S*randn(dim.p,1);
        end
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

% derive Free Energy at priors
F = FE(posterior,priors,suffStat,y);

if options.DisplayWin
    handles.hf = figure('name','VBA/MoG: convergence monitoring','color',[1 1 1]);
    for i=1:8
        handles.ha(i) = subplot(3,3,i,'parent',handles.hf,'box','off');
    end
    try; VBA_getSubplots ();end
    plot(handles.ha(1),F,'ro')
    set(handles.ha(1),'nextplot','add','ygrid','on')
    imagesc(sqrt(dist(y)),'parent',handles.ha(6))
    title(handles.ha(6),'distance data - data')
    set(handles.ha(7),'nextplot','add','ygrid','on')
    title(handles.ha(7),'data distribution')
    if dim.p >= 2
        [u,s,v,] = svd(y);
        yp = s(1:2,:)*v';
        plot(yp(1,:),yp(2,:),'.','color',0.8*[1 1 1],'parent',handles.ha(7))
        mup = u(:,1:2)'*posterior.muEta;
        handles.hp = plot(mup(1,:),mup(2,:),'k+','parent',handles.ha(7));
    else
        [ny,nx] = hist(y);
        bar(nx,ny./sum(ny),'parent',handles.ha(7),'facecolor',0.8*[1 1 1]);
        bounds = [min(y),max(y)];
        dy = diff(bounds)*1e-3;
        gridy = bounds(1):dy:bounds(2);
        pk = zeros(dim.K,length(gridy));
        Egam = EV(posterior);
        Er = posterior.d./sum(posterior.d);
        col = getColors(dim.K);
        for k=1:dim.K
            tmp = Egam(k)*(gridy-posterior.muEta(:,k)).^2;
            pk(k,:) = Er(k)*sqrt(2*pi*Egam(k)).*exp(-tmp/2);
            handles.hp(k) = plot(gridy,pk(k,:),'color',col(k,:),'parent',handles.ha(7));
        end
        set(handles.ha(7),'xlim',bounds)
    end
    VBA_getSubplots ();
end

% Main VB scheme
stop = 0;
it = 1;
while ~stop
    
    
    if options.DisplayWin % display results
        plot(handles.ha(1),it,F(it),'ro')
        title(handles.ha(1),'convergence monitoring')
        imagesc(posterior.z,'parent',handles.ha(2))
        title(handles.ha(2),'class labels')
        imagesc(sqrt(suffStat.Ed2),'parent',handles.ha(3))
        title(handles.ha(3),'distance data - components')
        imagesc(posterior.muEta,'parent',handles.ha(4))
        title(handles.ha(4),'components'' modes')
        cla(handles.ha(5))
        try
            plotUncertainTimeSeries(Egam,Egam./posterior.b_gamma,[],handles.ha(5));
        catch
            bar(Egam,'parent',handles.ha(5),'facecolor',0.8*[1 1 1])
        end
        title(handles.ha(5),'components'' precisions')
        set(handles.ha(5),'ygrid','on','box','off','xlim',[0.2,dim.K+0.8])
        if dim.p >= 2
            mup = u(:,1:2)'*posterior.muEta;
            set(handles.hp,'xdata',mup(1,:),'ydata',mup(2,:))
        else
            bounds = [min(y),max(y)];
            dy = diff(bounds)*1e-3;
            gridy = bounds(1):dy:bounds(2);
            pk = zeros(dim.K,length(gridy));
            Egam = EV(posterior);
            Er = posterior.d./sum(posterior.d);
            col = getColors(dim.K);
            for k=1:dim.K
                tmp = Egam(k)*(gridy-posterior.muEta(:,k)).^2;
                pk(k,:) = Er(k)*sqrt(2*pi*Egam(k)).*exp(-tmp/2);
                set(handles.hp(k),'ydata',pk(k,:),'color',col(k,:));
            end
        end
        [Ef,Vf] = VBA_dirichlet_moments(posterior.d);
        cla(handles.ha(8))
        try
            plotUncertainTimeSeries(Ef,Vf,[],handles.ha(8));
        catch
            bar(Ef,'parent',handles.ha(8),'facecolor',0.8*[1 1 1])
        end
        title(handles.ha(8),'components'' frequencies')
        set(handles.ha(8),'ygrid','on','box','off','xlim',[0.2,dim.K+0.8])
        drawnow
    end
    
    % ARD: components' death
    sumZ = posterior.z*ones(dim.n,1);
    k0 = find(sumZ < options.minSumZ);
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
        if dim.p < 2
            delete(handles.hp(k0))
            handles.hp(k0) = [];
        end
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
        posterior.SigmaEta{k} = VBA_inv(iS);
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
    
    % Convergence monitoring
    it = it +1;
    stop = checkStop(it,F,options);
    
end

% wrap up
if options.normalize
    out.normalize.Q = diag(sy.^2);
    out.normalize.m = my;
    % muEta <-- diag(sqrt(diag(Q)))*muEta + repmat(m,1,dim.K);
    % b_gamma(i) <-- = b_gamma(i).*Q;
end
out.dt = toc(options.tStart);
out.dim = dim;
out.options = options;
out.suffStat = suffStat;
out.date = clock;
out.F = F;
out.it = it;
out.F0 = LMEH0(y,out);
out.bor = 1/(1+exp(F(end)-out.F0));

if options.DisplayWin && dim.p >=2
    handles.ha(9) = subplot(3,3,9,'parent',handles.hf,'box','off');
    title(handles.ha(9),'2D-eigendata')
    out.handles = handles;
    VBA_projectMoG(posterior,out,y0);
    out = rmfield(out,'handles');
end

if options.verbose
    fprintf('---')
    fprintf('\n')
    fprintf(['Date: ',datestr(out.date),'\n'])
    if floor(out.dt./60) == 0
        timeString = [num2str(floor(out.dt)),' sec'];
    else
        timeString = [num2str(floor(out.dt./60)),' min'];
    end
    fprintf(['VB converged in ',num2str(it),' iterations (took ~',timeString,').','\n'])
    fprintf(['Dimensions:','\n'])
    fprintf(['     - observations: n=',num2str(dim.n),'\n'])
    fprintf(['     - features: p=',num2str(dim.p),'\n'])
    fprintf(['     - modes (final number): K=',num2str(dim.K),'\n'])
    fprintf(['Posterior probabilities:','\n'])
    fprintf(['     - MoG: p(H1|y)= ','%4.3f','\n'],1-out.bor)
    fprintf(['     - null: p(H0|y)= ','%4.3f','\n'],out.bor)
    fprintf('---')
    fprintf('\n')
    if options.normalize
        fprintf(['Warning: data has been normalized! (check out.normalize)','\n'])
        fprintf('---')
        fprintf('\n')
    end
end



% SUBFUNCTIONS

function stop = checkStop(it,F,options)
% checks stopping criteria
stop = 0;
if it<options.MinIter
    return
end
dF = F(it) - F(it-1);
if VBA_isWeird (F) || abs(dF)<=options.TolFun || it>=options.MaxIter
    stop = 1;
end

function F = FE(posterior,priors,suffStat,y)
% calculates Free Energy from sufficient statistics
[p,n] = size(y);
K = size(priors.muEta,2);
Elgam = ElogV(posterior);
Egam = EV(posterior);
Elogr = ElogR(posterior);
F = - n*p*log(2*pi)/2 - p*K*log(2*pi)/2 + p*K*log(2*pi*exp(1))/2 ...
    + gammaln(sum(priors.d)) - sum(gammaln(priors.d)) + (priors.d-1)'*Elogr ...
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
% calculates <log gamma>
Elv = psi(p.a_gamma) - log(p.b_gamma);

function Ev = EV(p)
% calculates <gamma>
Ev = p.a_gamma./p.b_gamma;

function Elr = ElogR(p)
% calculates <log r>
Elr = psi(p.d) - psi(sum(p.d));

function d2 = Ed2(y,p)
% calculates <(y_i-eta_k)'*(y_i-eta_k)>
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
    [~,i1,i2] = VBA_maxMat(d-diag(Inf*ones(size(eta,2),1)));
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

function d = dist(y)
n = size(y,2);
d = zeros(n,n);
for i=1:n
    tmp = repmat(y(:,i),1,n) - y;
    tmp = tmp.^2;
    d(i,:) = sum(tmp,1);
end

function colors = getColors(n)
hf = figure('visible','off');
ha = axes('parent',hf);
colors = get(ha,'colororder');
colors = repmat(colors,ceil(n/size(colors,1)),1);
colors = colors(1:n,:);
delete(hf)


function [F0] = LMEH0(y,out)
if out.dim.K == 1
    F0 = out.F(end);
else
    options.verbose = 0;
    options.DisplayWin = 0;
    [posterior,out] = VBA_MoG(y,1,options);
    F0 = out.F(end);
end
    


