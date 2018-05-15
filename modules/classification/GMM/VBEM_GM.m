function [xi,eta,F,out,K_opt] = VBEM_GM(y,K,options)

% function [xi,eta,F,theta,K_opt] = VBEM_GM(y,K,options)
% This function is a VB scheme for the MoG classifier.
% IN:
%   - y : the pxn data matrix, where n is the number of
%   'observations' to be classified, and p is the dimension of each of
%   these observations
%   - K : the maximum number of classes
%   - options : a structure containing infos regarding initialization and
%   inversion scheme (can be left empty -> defaults):
%       .TolFun = min increment of the free energy {1e-3}
%       .MaxIter = maximum number of iterations {200}
%       .MinIter = minimum number of iterations {50}
%       .verbose: flag for plotting the clustering results onto the
%       eigenspace of the data
%       .ARD : binary entry that enforces or not the component's death
%       process ({1})
%       .b : Dirichlet parameters for the class frequencies prior
%       .V : Variances of the components
%       .etaHat : Mean of the components
%       .Z : labelling process
%
% OUT:
%   - xi : the nx(K+1) matrix containing the estimated probability, for
%   each voxel, to belong to each of the clusters.
%   - eta : the (K+1)x3 matrix containing the estimated mean position of
%   the clusters
%   - F : free energy value series during the iteration process
%   - out : a structure variable containing the estimates of the
%   hyperparameters of the model, as well as posterior covariance
%   matrices and log evidence approximations.
%   - K_opt : Optimal dimension of the MoG, at convergence of the algorithm
%   (could be different of K because of the ARD component's death process).



%------------------ Dimensions of the model ----------------%
n = size(y,2); % number of observations to be classifed
p = size(y,1); % dimensions of the observations


%------------------ Maximization defaults -------------------------%
if ~exist('options','var') || ~isfield(options,'TolFun')
    options.TolFun = 1e-3;
end
if ~exist('options','var') || ~isfield(options,'MaxIter')
    options.MaxIter = 200;
end
if ~exist('options','var') || ~isfield(options,'MinIter')
    options.MinIter = 50;
end
if ~exist('options','var') || ~isfield(options,'verbose')
    options.verbose = 0;
end

%----------------- Initialization of the parameters -----------------%

% Sufficient statistics
Ytilde = y'*y;
% Data mean
mu = mean(y,2);

% Parameters of the Dirichlet prior pdf for the class frequencies
if exist('options','var') && isfield(options,'b')
    b = options.b;
else
    b = (1./(K)).*ones(K,1);
end
b00 = sum(b);
logZ0 = - gammaln(b)'*ones(K,1) + gammaln(b00);

% Component precisions
if exist('options','var') && isfield(options,'V')
    V = options.V;
else
    delta = K./var(y(:));
    V = (1/delta).*ones(K,1);
end

% Component means
if exist('options','var') && isfield(options,'etaHat')
    etaHat = options.etaHat;
else
    [etaHat] = dummyHierarchical(y,K);
%     for k = 1:K
%         etaHat(:,k) = mu + sqrt(V(k))*randn(p,1);
%     end
end

% Sufficient statistics of the components
Ntilde = diag(diag(etaHat'*etaHat) + p.*V);
T = y'*etaHat;

% Labelling process
if exist('options','var') && isfield(options,'Z')
    Z = options.Z;
else
    for k = 1:K
        lnPi(:,k) = psi(b(k)) - psi(b00) ...
            - (p./2).*log(V(k)) - (1./(2*V(k))).*( diag(Ytilde) + Ntilde(k,k) - 2.*T(:,k) );
    end
    [maxLnPi,indLnPi] = max(lnPi,[],2);
    Pi = exp(lnPi - repmat(maxLnPi,1,K));
    normFactor = sum(Pi,2);
    Z = Pi./repmat(normFactor,1,K);
end



%------------------------------------------------------------%
%----------------- Main VB learning scheme ------------------%
%------------------------------------------------------------%

stop = 0;
it = 1;
plotf = figure;
set(plotf,'menubar','none','numbertitle',...
    'off','name','Convergence monitoring')

while ~stop
    
    
    % counts for each class
    sumZ = Z'*ones(n,1);

    %----------------- ARD: component's death --------------%
    if any(sumZ == 0) && it > 1
        
        k0 = find(sumZ == 0);
        sumZ(k0) = [];
        K = length(sumZ);
        a_gamma(k0) = [];
        b_gamma(k0) = [];
        gammaHat(k0) = [];
        b(k0) = [];
        c(k0) = [];
        lambdaHat(k0) = [];
        V(k0) = [];
        etaHat(:,k0) = [];
        Ntilde = diag(diag(etaHat'*etaHat) + p.*V);
        T = y'*etaHat;
        Z(:,k0) = [];
        lnPi(:,k0) = [];
        b00 = sum(b);
        logZ0 = gammaln(b00) - gammaln(b)'*ones(K,1);
        
        disp([num2str(length(k0)),' component''s death: K_opt = ',num2str(K),'.'])
        
    end

    %---------------- MoG precision parameters ------------------%
    b_gamma = (1./2).*(Z'*diag(Ytilde) + Ntilde*sumZ - 2*diag(Z'*T));
    a_gamma = (p./2).*sumZ;
    gammaHat = a_gamma./b_gamma;
    
    %----- Class frequencies-----%
    c = b + sumZ;
    c0 = sum(c);
    lambdaHat = c./c0;

    %---- Mean of the clusters ----%
    previousEta = etaHat;
    for k = 1:K
        V(k) = 1./(gammaHat(k).*sumZ(k));
    end
    etaHat = y*Z;
    %     etaHat = etaHat*diag(gammaHat.*V);
    for k = 1:K
        etaHat(:,k) = etaHat(:,k)./sumZ(k);
    end

    Ntilde = diag(diag(etaHat'*etaHat) + p.*V);
    T = y'*etaHat;

    %----- Labelling process -----%
    for k = 1:K
        lnPi(:,k) = psi(c(k)) - psi(c0) ...
            + (p./2).*(psi(a_gamma(k)) - log(b_gamma(k))) ...
            - 0.5.*gammaHat(k).*( diag(Ytilde) + Ntilde(k,k) - 2.*T(:,k) );
    end
    [maxLnPi,indLnPi] = max(lnPi,[],2);
    Pi = exp(lnPi - repmat(maxLnPi,1,K));
    normFactor = sum(Pi,2);
    Z = Pi./repmat(normFactor,1,K);

    %------- Calculation of the log-evidence ------%
    logZ = gammaln(c)'*ones(K,1) - gammaln(c0);
    i0 = find(Z==0);
    ZZ = Z;
    ZZ(i0) = [];
    
    F(it) = sum(sum(Z.*lnPi)) - (p/2)*log(2*pi)*sum(sum(Z)) ... % <ln p(y,xi|mu,lambda,gamma)>
        + sum((b-1).*(psi(c) - psi(c0))) + logZ0 ...            % <ln p(lambda)>
        - sum(psi(a_gamma) - log(b_gamma)) ...                  % <ln p(gamma)>
        - sum(sum(ZZ.*log(ZZ))) ...                             % H(q(xi))
        + (K*p/2).*(log(2*pi)+1) + (p/2)*log(V)'*ones(K,1) ...  % H(q(mu))
        + ( gammaln(a_gamma) - log(b_gamma) + ...
        (1-a_gamma).*psi(a_gamma) + a_gamma )'*ones(K,1) ...    % H(q(gamma))
        + logZ - ((c-1).*(psi(c)-psi(c0)))'*ones(K,1);          % H(q(lambda))

    %------ Graphical output and convergence monitoring -----%
    figure(plotf);
    subplot(2,1,1),hold on,plot(it,F(it),'ro'),drawnow
    
    if it > 1
        dF(it) = F(end)-F(end-1);
        figure(plotf);
        subplot(2,1,2),hold on,plot(it,dF(it),'ro'),drawnow
        if it >= options.MinIter  && ...
                (abs(dF(end)) < options.TolFun ||it >= options.MaxIter)
            stop  = 1;
        end
    else
        figure(plotf);
        subplot(2,1,1),title('Variational free energy (F)')
        subplot(2,1,2),title('First derivative of F')
    end
    
    it = it +1;
    
end



%--------- Output formulation ---------%

xi = Z;
eta = etaHat;
K_opt = K;
out.gamma = gammaHat;
out.varEta = V;
out.lambda = lambdaHat;
out.varLambda = c.*(c0-c)./(c0.^2.*(c0+1));

%------ Approximations to the model evidence --------%

% Free energy
out.logEv.F = F(end);

% Fuzzy hyper-volume
out.logEv.logFHV = -(p/2).*log(gammaHat)'*ones(K,1);

% Evidence density
for k = 1:K
    for i = 1:n
        L(i,k) = -(p/2)*log(2*pi) ...
            + (p./2).*log(gammaHat(k)) ...
            - (1/2).*gammaHat(k).*sum((y(:,i)-etaHat(:,k)).^2);
    end
end
P = L.*Z;
LogLikelihood = sum(sum(P));
out.logEv.logRho = LogLikelihood./out.logEv.logFHV;

% Minimum description length
out.logEv.MDL = - LogLikelihood + (1/2).*(n+2+p)*K*log(n);

% Partition coefficient
P_lambda = logZ0 + (b-1)'*log(lambdaHat);
PP = P + P_lambda;
out.logEv.PC = (1/n).*sum(sum(exp(2*PP)));


% Restricted likelihood
i0 = find(lnPi==0);
lnPi(i0) = [];
Z(i0) = [];
out.logEv.ReL = sum(sum(Z.*lnPi))  - (n+2+p)*K*log(2*pi);



%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%




function [f] = psi(z)

siz = size(z);
z=z(:);
zz=z;
f = 0.*z; % reserve space in advance
%reflection point
p=find(real(z)<0.5);
if ~isempty(p)
    z(p)=1-z(p);
end
%Lanczos approximation for the complex plane
g=607/128; % best results when 4<=g<=5

c = [  0.99999999999999709182;
    57.156235665862923517;
    -59.597960355475491248;
    14.136097974741747174;
    -0.49191381609762019978;
    .33994649984811888699e-4;
    .46523628927048575665e-4;
    -.98374475304879564677e-4;
    .15808870322491248884e-3;
    -.21026444172410488319e-3;
    .21743961811521264320e-3;
    -.16431810653676389022e-3;
    .84418223983852743293e-4;
    -.26190838401581408670e-4;
    .36899182659531622704e-5];
n=0;
d=0;
for k=size(c,1):-1:2
    dz=1./(z+k-2);
    dd=c(k).*dz;
    d=d+dd;
    n=n-dd.*dz;
end
d=d+c(1);
gg=z+g-0.5;
%log is accurate to about 13 digits...
f = log(gg) + (n./d - g./gg) ;
if ~isempty(p)
    f(p) = f(p)-pi*cot(pi*zz(p));
end
p=find(round(zz)==zz & real(zz)<=0 & imag(zz)==0);
if ~isempty(p)
    f(p) = Inf;
end
f=reshape(f,siz);
return

function [eta,Z] = dummyHierarchical(y,K)
[p,n] = size(y);
Z = eye(n);
eta = y;
while size(Z,1) > K
    C = corrcoef(eta);
    [~,i,j] = VBA_maxMat(C-diag(Inf*ones(n,1)));
    newZ = eye(n);
    newZ(i,j) = 1;
    newZ(j,:) = [];
    Z = newZ*Z;
    eta = y*Z'*diag(sum(Z,2).^-1);
    n = size(eta,2);
end

