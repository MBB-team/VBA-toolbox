function [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity(f_fname,g_fname,u,n_t,options,dim,N,np,lx,ly)

% prior predictive density under sDCM generative model (MCMC)
% function [pX,gX,pY,gY,X,Y] =
% VBA_MCMC_predictiveDensity(f_fname,g_fname,u,n_t,options,dim,N,np)
% IN:
%   - f_fname: name/handle of the evolution function
%   - g_fname: name/handle of the observation function
%   - u: input to the system
%   - n_t: number of time samples of the predictive density
%   - options: the optional structure. The prior pdfs on
%   evolution/observation parameters/hyperparameters are extracted from
%   this structure.
%   - dim: the dimension of the sDCM generative model
%   - N: the number of MCMC samples used in the derivatio of the predictive
%   density
%   - np: the number of bins used to form the MCMC empirical histograms
%   - lx: dim.nX2 matrix of min/max admissible range of states values
%   - ly: dim.pX2 matrix of min/max admissible range of observables values
% OUT:
%   - pX: the n_tXnpXdim.n 3D array of MCMC empirical histograms of hidden
%   states
%   - gX: the npXdim.n 2D array giving the grid used for forming the MCMC
%   empirical histograms on each dimension of the hidden states
%   - pX/gY: [id, but for observed data]
%   - X: the dim.nXn_tXN set of MCMC samples used to form the histograms
%   - Y: [id, but for observed data]


% Get time
et0 = clock;

% default sample size and histogram resolution
try, N ; catch, N  = 1e3 ; end
try, np; catch, np = 50  ; end
try, lx; catch, lx = []  ; end
try, ly; catch, ly = []  ; end


% fix precision parameters
alpha = options.priors.a_alpha./options.priors.b_alpha;
sigma = options.priors.a_sigma./options.priors.b_sigma;

% fill in missing optional fields
try
    if isinf(options.priors.a_alpha)
        options.priors.a_alpha = 1;
        options.priors.b_alpha = 1;
    end
end
options.verbose = 0;
dim.n_t = n_t;
[options] = VBA_check([],u,f_fname,g_fname,dim,options);

fprintf(1,'MCMC sampling...')
fprintf(1,'%6.2f %%',0)
Y = zeros(dim.p,n_t,N);
X = zeros(dim.n,n_t,N);
out = [];
for i=1:N
    [x0,theta,phi] = sampleFromPriors(options,dim);
    try
        [y,x] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);
        ok = ~ VBA_isWeird ({x, y}) && VBA_isInRange (x, lx) && VBA_isInRange (y, ly);
    catch
        ok = false;
    end
    if ok
        Y(:,:,i) = y;
        X(:,:,i) = x;
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',i*100/N)
    else
        out = [out;i];
    end
end
Y(:,:,out) = [];
X(:,:,out) = [];
fprintf(1,repmat('\b',1,8))
fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
fprintf(1,'\n')

et0 = clock;
fprintf(1,'Forming histograms along X dimensions...')
fprintf(1,'%6.2f %%',0)
pX = zeros(n_t,np,dim.n);
gX = zeros(np,dim.n);
for i=1:dim.n
    Xi = squeeze(X(i,:,:));
    if isempty(lx)
        m = mean(Xi(:));
        sv = std(Xi(:));
        nx = m-3*sv:6*sv/(np-1):m+3*sv;
    else
        dx = lx(i,2)-lx(i,1);
        nx = lx(i,1):dx/(np-1):lx(i,2);
    end
    [ny,nx] = hist(Xi',nx);
    pX(:,:,i) = ny';
    gX(:,i) = nx;
    fprintf(1,repmat('\b',1,8))
    fprintf(1,'%6.2f %%',i*100/dim.n)
end
fprintf(1,repmat('\b',1,8))
fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
fprintf(1,'\n')

et0 = clock;
fprintf(1,'Forming histograms along Y dimensions...')
fprintf(1,'%6.2f %%',0)
pY = zeros(n_t,np,dim.p);
gY = zeros(np,dim.p);
for i=1:dim.p
    Yi = squeeze(Y(i,:,:));
    if isempty(ly)
        m = mean(Yi(:));
        sv = std(Yi(:));
        nx = m-3*sv:6*sv/(np-1):m+3*sv;
    else
        dy = ly(i,2)-ly(i,1);
        nx = ly(i,1):dy/(np-1):ly(i,2);
    end
    [ny,nx] = hist(Yi',nx);
    pY(:,:,i) = ny';
    gY(:,i) = nx;
    fprintf(1,repmat('\b',1,8))
    fprintf(1,'%6.2f %%',i*100/dim.p)
end
fprintf(1,repmat('\b',1,8))
fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
fprintf(1,'\n')



function [x0,theta,phi] = sampleFromPriors(options,dim)

priors = options.priors;

if dim.n > 0
    if ~isequal(priors.SigmaX0,zeros(size(priors.SigmaX0)))
        sV = VBA_sqrtm (priors.SigmaX0);
        x0 = priors.muX0 + sV*randn(dim.n,1);
    else
        x0 = priors.muX0;
    end
else
    x0 = [];
end

if dim.n_theta > 0
    if ~isequal(priors.SigmaTheta,zeros(size(priors.SigmaTheta)))
        sV = VBA_sqrtm(priors.SigmaTheta);
        theta = priors.muTheta + sV*randn(dim.n_theta,1);
    else
        theta = priors.muTheta;
    end
else
    theta = [];
end

if dim.n_phi > 0
    if ~isequal(priors.SigmaPhi,zeros(size(priors.SigmaPhi)))
        sV = VBA_sqrtm (priors.SigmaPhi);
        phi = priors.muPhi + sV*randn(dim.n_phi,1);
    else
        phi = priors.muPhi;
    end
else
    phi = [];
end


