function demo_sparsityPrior ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% demo_sparsityPrior ()
% demo of penalized regression (Bayesian LASSO)
%
% This demo shows how to perform a GLM regression using a L1-norm, L2-norm, 
% or an adaptive norm regularization.
% Critically, we vary the number of "dummy" regressors (that are not used
% to generate the data), and assess how this level of sparsity affects the
% accuracy of the respective regularized regressions.
%
% Background:
% ~~~~~~~~~~~
%
% When the number of parameters to estimate is higher than the number of
% data points, it can be usefull to apply sparsity constraints during model
% fitting. In general, this is done by defining a penalty term (eg the l1-
% or l2-norm) that help shrinking the number of effective parameters, as
% for example in the LASSO estimator or elastic nets. 
% Here, we demonstrate how applying a simple transformation on the
% parameters allows us to emulate l1-regularization in a regression, and
% therefore perform a Bayesian LASSO estimation.
%
% See arXiv:1703.07168 for more details.
%
% /////////////////////////////////////////////////////////////////////////

% Illustrating the transformation
% =========================================================================
displaySparsifyTransform ()

% Parameters of the demo
% =========================================================================
% number of observations (data points)
N = 16; 
% total number of regressors
K = 32;
% number of useless regressors (0 weight)
sparsity = 16 : 3 : K ; 
% number of Monte - Carlo simulations for each sparsity level
M = 30; 

% observation precision
sigma = 10;

% options
options.DisplayWin = 0;
options.verbose = 0;

% main loop
% =========================================================================
% loop over different levels of sparsity underlying the data
for i = 1 : numel(sparsity)
    % do multiple times
    for m = 1 : M
        
        % simulate data
        % -----------------------------------------------------------------
        % design matrix (set or regressors)
        X = randn (N, K);
        % weight of the regressors
        phi = [ones(K - sparsity(i), 1) ; zeros(sparsity(i), 1)]; 
        % generate random observations
        y = X * phi +  VBA_random ('Gaussian', 0,  1 / sigma, N, 1);
        % store for inversion
        options.inG.X = X;
        
        % L1 estimator (emulated Laplace priors)
        % -----------------------------------------------------------------
        dims.n_phi = K;
        [p1, o1] = VBA_NLStateSpaceModel (y, [], [], @g_GLMsparse, dims, options);
        % recover transformed parameters
        phiHat = VBA_sparsifyPrior (p1.muPhi);
        accuracy(2, i, m) = myCorr (phi, phiHat);
        
        % L2 estimator (usual Gaussian priors)
        % -----------------------------------------------------------------
        dims.n_phi = K;
        [p2, o2] = VBA_NLStateSpaceModel (y, [], [], @g_GLM, dims, options);
        phiHat = p2.muPhi;
        accuracy(1, i, m) = myCorr (phi, phiHat);
        
        % adaptive sparse-estimator
        % -----------------------------------------------------------------
        % one aditional parameter for the exponent
        dims.n_phi = K + 1;
        % priors set to l1-norm
        optionsAdapt = options;
        optionsAdapt.priors.muPhi = [zeros(K, 1); log(2)];
        optionsAdapt.priors.SigmaPhi = diag([ones(K, 1); 0.2]);
        % estimation
        [p3, o3] = VBA_NLStateSpaceModel (y, [], [], @g_GLMsparseAdapt, dims, optionsAdapt);
        % get estimated exponent
        exponent(i, m) =  exp (p3.muPhi(K + 1));
        % recover transformed parameters
        phiHat = VBA_sparsifyPrior (p3.muPhi(1 : end - 1), log (exponent(i, m)));
        accuracy(3, i, m) = myCorr (phi, phiHat);

        
    end
end

% Display
% =========================================================================
hf = VBA_figure ('name', 'demo_sparsityPrior: results');

% accuracy
ha = subplot(1,2,1,'parent',hf);
mr = mean(accuracy,3);
sr = std(accuracy,[],3) / sqrt (M);
hb = errorbar(ha,mr',sr');
set(hb(1), 'Color', 'b')
set(hb(2), 'Color', 'r')
set(hb(3), 'Color', 'm')
legend(ha,{'L2-norm (Gaussian)','L1-norm (Laplace)','Adaptive sparsity'},'Location','northwest')
xlabel(ha,'sparsity of generative model')
ylabel(ha,'parameter estimation accuracy')
box off
% adaptive exponent
ha = subplot(1,2,2,'parent',hf);
mse = mean(exponent,2);
sse = std(exponent,[],2) / sqrt (M);
hb = errorbar(ha,mse,sse);
set(hb, 'Color', 'm');
xlabel(ha,'sparsity of generative model')
ylabel(ha,'estimated sparsity exponent')
box off
end

% Subfunctions
% #########################################################################
function displaySparsifyTransform ()

hf = VBA_figure ('name', 'demo_sparsityPrior: sparsify mapping');

x = -4 : 0.01 : 4;
gridP = [- log(2), 0, log(2)];
for i = 1 : length (gridP)
    sx(i, :) = VBA_sparsifyPrior (x,  gridP(i));
end
ha = subplot (1, 2, 1, 'parent', hf);
hb = plot (ha, x, sx');
set(hb(1), 'Color', 'g')
set(hb(2), 'Color', 'b')
set(hb(3), 'Color', 'r')
box off
xlabel ('original parameter')
ylabel ('transformed parameter')
legend (ha, ...
    {'sparsity exponent = 1/2', 'sparsity exponent = 1 ( l2-norm )', 'sparsity exponent = 2 ( l1-norm )' }, ...
    'Location', 'northwest');


ha = subplot (1, 2, 2, 'parent', hf);
for i = 1 : length (gridP)
    ps = interp1 (sx(i, :), normpdf (x, 0, 1), x);
    ps = ps / nansum (ps);
    hb(i) = plot (ha, x, ps);
    hold on;
end
set(hb(1), 'Color', 'g')
set(hb(2), 'Color', 'b')
set(hb(3), 'Color', 'r')
hold off
box off
xlabel ('transformed parameter')
ylabel ('prior density')



end

function r = myCorr (x, y)
    tmp = corrcoef (x, y);
    r = tmp(1, 2);
end

