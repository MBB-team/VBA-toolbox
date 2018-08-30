function [post,out] = VB_CRP(y,alpha,options)
% function [post,out] = VB_CRP(y,alpha,options)
% VB scheme for the CRP-MoG classifier.
% IN:
%   - y : the pxn data matrix, where n is the number of
%   'observations' to be classified, and p is the dimension of each of
%   these observations
%   - alpha : CRP concentration parameter
%   - options : a structure containing infos regarding initialization and
%   inversion scheme (can be left empty -> defaults):
%       .priors : structure of priors for generative density
%           .mu0: px1 prior expectation of Gaussian component mean
%           .S0: prior covariance matrix of Gaussian component mean
%           .a0: prior counts for each class
%           .b0: prior sum of squared error within each class
%           NB: a0./b0 is the a priori expected precision (inverse
%           variance) of Gaussian components.
%       .TolFun = min increment of the free energy {1e-3}
%       .MaxIter = maximum number of iterations {128}
%       .MinIter = minimum number of iterations {1}
%       .verbose = flag for verbose mode {0}
%
% OUT:
%   - post : posterior structure:
%       .xi: Kxn matrix of labels probabilities (K=estimated number of
%       classes)
%       .mu: pxK matrix of component estimated mean
%       .S: Kx1 cell array of posterior covariance matrices over component
%       means
%       .a: Kx1 vector of estimated counts for each class
%       .b: Kx1 vector of etimated sum of squared error within each class
%       NB: a./b gives the expected precision (inverse variance) of each
%       Gaussian component!
%   - out : output structure containing:
%       .xih: the series of estimated labels
%       .F: lower bound on the log model evidence

disp('-- Variational Bayesian CRP classification --')

[p,n] = size(y); % dimensions of the observations

% options structure
if ~exist('alpha','var') || isempty(alpha)
    alpha = 1;
end
if ~exist('options','var') || ~isfield(options,'priors')
    options.priors = struct(...
        'mu0',zeros(p,1),...
        'S0',eye(p),...
        'a0',1,...
        'b0',1,...
        'alpha',alpha);
end
if ~exist('options','var') || ~isfield(options,'TolFun')
    options.TolFun = 1e-2;
end
if ~exist('options','var') || ~isfield(options,'MaxIter')
    options.MaxIter = 64;
end
if ~exist('options','var') || ~isfield(options,'MinIter')
    options.MinIter = 1;
end
if ~exist('options','var') || ~isfield(options,'verbose')
    options.verbose = 0;
end


% pre-allocate sufficient statistics
post.xi = zeros(1,n);
post.mu = [];
post.S = cell(0);
post.a = [];
post.b = [];
iS0 = pinv(options.priors.S0);
m0 = options.priors.mu0;
iS0m0 = iS0*m0;
Ip = eye(p);
yy = diag(y'*y);

% Initial condition
post.xi(1) = 1;
post.a(1) = options.priors.a0;
post.b(1) = options.priors.b0;
stop = 0;
if options.verbose
    hf = figure('color',[1 1 1]);
    ha(1) = subplot(2,1,1,...
        'parent',hf,...
        'nextplot','add');
    ha(2) = subplot(2,1,2,...
        'parent',hf,...
        'nextplot','add');
    drawnow
end
F = zeros(1,options.MaxIter);
F(1) = -Inf;
it = 1;
while ~stop
    % approximate posterior over 1stOM of MoG components
    post.S{1} = pinv(...
        iS0 + (post.a(1)./post.b(1))*Ip*post.xi(1)...
        );
    post.mu(:,1) = post.S{1}*(...
        iS0m0 + (post.a(1)./post.b(1))*post.xi(1)*y(:,1)...
        );
    % approximate posterior over 2ndOM of MoG components
    dy = y(:,1) - post.mu(:,1);
    post.a(1) = options.priors.a0 + (p/2)*post.xi(1);
    post.b(2) = options.priors.b0 + ...
        (1/2)*post.xi(1)*(trace(dy'*dy)+trace(post.S{1}));
    
    % convergence
    it = it + 1;
    F(it) = getFE(1,y,post,options.priors);
    dF = F(it) - F(it-1);
    if dF <= options.TolFun
        stop = 1;
    end
    if options.verbose
        plot(ha(1),it,F(it),'r+')
        plot(ha(2),it,dF,'r+')
        axis(ha,'tight')
        drawnow
    end
    
end

% Initialize accumulation variables
counts = sum(post.xi,2);
Zy = y*post.xi';


% Loop over sequence of data points
fprintf(1,'Looping over samples ...')
fprintf(1,'%6.2f %%',0)

for i=2:n
    
    % prior predictive belief over label from CRP
    post.xi = [post.xi;zeros(1,n)];
    post.xi(:,i) = [counts;alpha]./(alpha+i-1);
    
    % remove components that are less probable than a new one
    empty = find(counts<alpha./(alpha+i-1));
    post = killCRP(post,empty);
    
    % initialize accumulation variables
    counts = sum(post.xi,2);
    Zy = y*post.xi';
    
    % initialize belief about new component...
    Ki = size(post.xi,1);
    post.a(Ki) = options.priors.a0;
    post.b(Ki) = options.priors.b0;
    post.mu(:,Ki) = m0;
    post.S{Ki} = options.priors.S0;
    
    % run VB update rules
    F = zeros(1,options.MaxIter);
    F(1) = -Inf;
    stop = 0;
    iti = 0;
    while ~stop
        lA = zeros(Ki,1);
        for k=1:Ki
            % approximate posterior over 1stOM of MoG components
            post.S{k} = pinv(iS0 + (post.a(k)./post.b(k))*Ip*counts(k));
            post.mu(:,k) = ...
                post.S{k}*( iS0m0 + (post.a(k)./post.b(k))*Zy(:,k) );
            % approximate posterior over 2ndOM of MoG components
            ymu = y'*post.mu(:,k);
            mu2 = post.mu(:,k)'*post.mu(:,k);
            dy = y(:,1:i) - repmat(post.mu(:,k),1,i);
            post.a(k) = options.priors.a0 + (p/2)*counts(k);
            post.b(k) = options.priors.b0...
                +(1/2)*post.xi(k,1:i)*(yy(1:i) - 2*ymu(1:i) + mu2)...
                +(1/2)*trace(post.S{k})*counts(k);
            % log-Gaussian responsibility
            lA(k) = ...
                -(p/2)*log(2*pi) ...
                -(1/2)*(post.a(k)./post.b(k))*(dy(:,i)'*dy(:,i)+trace(post.S{k}))...
                -(p/2)*(psi(post.a(k)))-log(post.b(k));
        end
        % approximate posterior over labels
        pp = lA + log(counts) - 0.5.*sum(post.xi.*(1-post.xi),2)./counts.^2;
        pp = exp(pp - max(pp));
        pp(isnan(pp)|isinf(pp)) = 0;
        post.xi(:,i) = pp./sum(pp);
        
        % update accumulation variables
        counts = sum(post.xi,2);
        Zy = y*post.xi';
        
        % convergence
        it = it + 1;
        iti = iti + 1;
        F(it) = getFE(1,y,post,options.priors);
        dF = F(it) - F(it-1);
        if hasconv(dF,iti,options)
            stop = 1;
        end
        
        % display
        if options.verbose
            plot(ha(1),it,F(it),'r+')
            plot(ha(2),it,dF,'r+')
            if stop
                plot(ha(1),it,F(it),'r*')
                plot(ha(2),it,dF,'r*')
            end
            axis(ha,'tight')
            drawnow
        end
        fprintf(1,'\b\b\b\b\b\b\b\b')
        fprintf(1,'%6.2f %%',100*i/n)
        
    end % VB update loop
    
end % loop over samples

fprintf(1,'\b\b\b\b\b\b\b\b')
fprintf(' OK.')
fprintf('\n')

% remove components that are less probable than a new one
empty = find(counts<alpha./(alpha+i-1));
post = killCRP(post,empty);

out.xih = zeros(1,n);
for i=1:n
    out.xih(i) = find(post.xi(:,i)==max(post.xi(:,i)));
end
out.F = F(end);




function Fi = getFE(i,y,post,priors)
[p,n] = size(y);
K = size(post.xi,1);
lA = zeros(K,1);
iS0 = pinv(priors.S0);
Seta = 0;
Ssig = 0;
Sz = 0;
Ep = 0;
Epz = 0;
for k=1:K
    Seta = Seta + entropyGauss(post.S{k});
    Ssig = Ssig + entropyGamma(post.a(k),post.b(k));
    dm = priors.mu0 - post.mu(:,k);
    Ep = Ep...
        + (priors.a0-1)*(psi(post.a(k)))-log(post.b(k))...
        - priors.b0.*post.a(k)./post.b(k)...
        - 0.5*(dm'*iS0*dm+trace(post.S{k}*iS0));
    dy = y(:,i) - post.mu(:,k);
    lA(k) = ...
        -(p/2)*log(2*pi) ...
        - 0.5*(...
        post.a(k)./post.b(k)*(dy'*dy+trace(post.S{k}))...
        -p*(psi(post.a(k)))-log(post.b(k)));
    if i>1
        counts = sum(post.xi(k,1:i-1),2);
        Epz = Epz ...
            + post.xi(k,i)*(...
            log(counts) - 0.5.*sum(post.xi.*(1-post.xi),2)./counts.^2)...
            - log(priors.alpha+i-1);
        Sz = Sz + entropyBinomial(post.xi(:,i));
    end
end
Fi = post.xi(:,i)'*lA + Epz + Ep + Seta + Ssig + Sz;

function S = entropyGauss(S)
n = size(S,1);
S = log(sqrt(det(S)*(2*pi*exp(1)).^n));

function S = entropyGamma(a,b)
S = a - log(b) + gammaln(a) + (1-a)*psi(a);

function S = entropyBinomial(p)
p(p==0) = 1e-3;
p(p==1) = 1-1e-3;
S = sum(p.*log(p));

function post = killCRP(post,ind)
if ~isempty(ind)
    post.xi(ind,:) = [];
    sp = sum(post.xi,1);
    sp(sp==0) = 1;
    post.xi = post.xi./repmat(sp,size(post.xi,1),1);
    post.a(ind) = [];
    post.b(ind) = [];
    post.mu(:,ind) = [];
    post.S{ind} = [];
end

function x = hasconv(dF,it,options)
x = 0;
if it<options.MinIter
    return
else
    if abs(dF) <= options.TolFun ||...
            it>=options.MaxIter
        x = 1;
    end
end


