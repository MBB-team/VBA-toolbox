function [all] = VBA_classification(X,y,k,verbose,options,sparse)
% performs binary classification using VBA
% function [pv,stat,df,all] = VBA_classification(X,y,k,verbose,options)
% In brief, this function fits the following logistic regresion model:
%   p(y=1) = s(X*beta)
% where y is the binary data, s is the standard sigmoid mapping, X is the
% design matrix and beta are unknown weight parameters.
% Statistical testing is performed using cross-validation scheme, such that
% the p-value computes the probability of the test score under the null H0.
% In cross-validation, the test score is the binary accuracy outcome
% sampled across test sets. The null distribution accounts for potential
% imbalance in the data in that the first-order moment of the corresponding
% binomial distribution is specified in terms of the sample mean of the
% data. VBA_classification uses a k-fold cross-validation scheme, which
% involves partitioning the data into k subsamples of equal size, each of
% which is used as the test set in turn.
% Note: VBA_classification also performs a fully Bayesian model inversion
% on the entire dataset (full-data inversion).
% IN:
%   - X: nXp design matrix
%   - y: nX1 binary data matrix
%   - k: the number of folds for the k-fold cross-validation scheme. If k
%   is set to n (default), then VBA_classification uses a leave-one-out
%   scheme. If k is set to 0, then the cross-validation scheme is entirely
%   skipped (only the full-data inversion is performed).
%   - verbose: flag for displaying results (default is 0)
%   - options: structure containing VBA's options (can be used for passing,
%   e.g., priors on classification weight parameters)
%   - sparse: when sparse=1, VBA_classification uses sparsifying priors
%   (default=0).
% OUT:
%   - all: structure array with fields:
%       .stat: structure of summary statistics:
%           .pv: classical p-value on classifier accuracy
%           .success: nx1 vector of successful classifications
%           .pa: cross-validation prediction accuracy
%           .bpa: balanced prediction accuracy
%           .pBayes: Bayesian exceedance prob. on prediction accuracy
%       .P: pxk matrix of estimated classification weights
%       .r: data imbalance (r=0.5 means balanced data)
%       .in: structure storing the inputs to VBA_classification
%       .date: date stamp
%       .dt: computing time
%       .handles: a structure containing the handles of the graphical
%       objects (filled-in by VBA_classication_display.m)
% Note: cross-validation results can be displayed using the following
% command line: VBA_classification_display(all).

% fill in default I/O
tStart = tic;
all = [];
[n,m] = size(y);
[n0,p] = size(X);

try,k;catch,k=n;end
try,verbose;catch,verbose=0;end
try % check whether VBA's options structure has been sepecified
    options;
catch
    options = [];
end
try % use sparsify tansform?
    sparse = ~~sparse(1);
catch
    sparse = 0;
end

% check basic numerical requirements
try
    if VBA_isWeird (y)
        disp('Error: data contains weird values!')
        return
    end
end
if ~ VBA_isBinary (y)
    disp('Error: data should be binary!')
    return
end
if ~isequal(n,n0)
    disp('Error: design matrix has to have as many rows as the data matrix!')
    return
end
if ~isequal(m,1)
    disp('Error: data should be a nx1 vector!')
    return
end
if k>n || k<0
    disp('Warning: number of folds reduced to data length!')
    k = n;
end


% randomly re-order the data (for k-fold partitioning)
if k ~= n
    io = randperm(n);
    y = y(io);
    X = X(io,:);
end

% specify options and priors for VBA inversion
dim.p = n;
dim.n_t = 1;
dim.n_phi = p;
dim.n_theta = 0;
dim.n = 0;
g_fname = @g_classif0;
options.sources = struct('type',1,'out',1);
options.DisplayWin = 0;
options.verbose = 0;
options.inG.X = X';
options.inG.sparse = sparse;
options.n0 = 0; % number of dummy counts

if ~isequal(k,0) % performing cross-validation scheme
    if verbose
        disp(' ')
        fprintf(1,'Performing k-folds cross-validation scheme...')
        fprintf(1,'%6.2f %%',0)
        et0 = clock; % get time
    end
    sizeFolds = floor(n./k);
    acc = zeros(n,1); % test accuracy
    acc0 = zeros(n,1); % training accuracy
    Eg = zeros(n,1); % out-of-sample prediction E[y] (on test data)
    Vg = zeros(n,1); % out-of-sample prediction V[y] (on test data)
    P = zeros(p,k);
    for i=1:k
        if i<k
            itest = (i-1)*sizeFolds+1:i*sizeFolds;
        else % last test fold contains all remaining data
            itest = (i-1)*sizeFolds+1:n;
        end
        options.isYout = zeros(n,1);
        options.isYout(itest) = 1;
        [posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
        Eg(itest) = out.suffStat.gx(itest);
        Vg(itest) = out.suffStat.vy(itest);
        ytest = Eg(itest)>=0.5;
        acc(itest) = [y(itest)==ytest];  
        acc0(itest) = out.fit.acc;
        if sparse
            P(:,i) = VBA_sparsifyPrior (posterior.muPhi);
        else
            P(:,i) = posterior.muPhi;
        end
        if verbose
            fprintf(1,repmat('\b',1,8))
            fprintf(1,'%6.2f %%',floor(100*i/k))
        end
    end
    if verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
        fprintf(1,'\n')
    end
end

% performing full-data inversion
if verbose
    fprintf(1,'Performing whole-data inversion...')
    et0 = clock; % get time
end
options.isYout = zeros(n,1);
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
if verbose
    fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
    fprintf(1,'\n')
end

% wrap-up
all.in.X = X;
all.in.y = y;
all.in.k = k;
all.in.sparse = sparse;
all.date = clock;
all.dt = toc(tStart);
all.posterior = posterior;
all.out = out;
if ~isequal(k,0)
    r = mean(y);
    if r<0.5
        r = 1-r;
    end
    nok = sum(acc);
    [pdf0,cdf0] = VBA_binomial(nok,n,r);
    all.stat.pv = 1 - cdf0; % classical p-value on classifier accuracy
    all.stat.success = acc; % nx1 vector of successful classifications
    all.stat.pa = nok./n; % cross-validation prediction accuracy
    all.stat.bpa = 0.5*(sum(acc.*y)./sum(y) + sum(acc.*(1-y))./sum(1-y)); % balanced class. acc.
    all.stat.pBayes = VBA_PPM(nok+options.n0,n-nok+options.n0,r,'beta',0); % Bayesian exceedance prob.
    all.Eg = Eg;
    all.Vg = Vg;
    all.P = P; % set of estimated weights (for each k-fold)
    all.r = r; % data imbalance (r=0.5 means balanced data)
    if verbose
        [all] = VBA_classification_display(all);
    end
end


