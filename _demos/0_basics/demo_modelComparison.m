function demo_modelComparison ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% demo_modelComparison ()
% demo of random-effect group-level bayesian model selection in the context
% of nested linear models. 
%
% /////////////////////////////////////////////////////////////////////////

% definition of the observations and models
% =========================================================================

% dimensions
% -------------------------------------------------------------------------
% number of subjects
nSubjects = 32;
% number of observatiosn per subject
nObs = 64;
% number of predictors (regressors)
nPred = 4;

% parameters
% -------------------------------------------------------------------------
% noise variance
sigma2 = 1e0;
% signal strength
signal = 1e0; 

% generate data
% -------------------------------------------------------------------------
for i = 1 : nSubjects
    
    % full design matrix
    X1 = randn (nObs, nPred);
    % nested model
    X2 = X1(:, 1 : 2);
    
    % random weights
    b = sqrt (signal) * rand(nPred, 1);
    
    % simulate observations
    y1 = X1 * b + sqrt (sigma2) * randn (nObs, 1);
    y2 = X2 * b(1 : 2) + sqrt (sigma2) * randn (nObs, 1);
    
    % compute model evidence (frequentist limit)
    logEvidence_y1(1, i) = lev_GLM (y1, X1);
    logEvidence_y1(2, i) = lev_GLM (y1, X2);
    
    logEvidence_y2(1, i) = lev_GLM (y2, X1);
    logEvidence_y2(2, i) = lev_GLM (y2, X2);
    
end

% display empirical histogram of log-Bayes factors
% -------------------------------------------------------------------------
plotBayesFactor (logEvidence_y1, logEvidence_y2);

% perform model selection with the VBA
% =========================================================================
options.verbose = false;

% perform group-BMS on data generated under the full model
[p1, o1] = VBA_groupBMC (logEvidence_y1, options);
set (o1.options.handles.hf, 'name', 'group BMS: y_1')

fprintf('Statistics in favor of the true model (m1): pxp = %04.3f (Ef = %04.3f)\n', o1.pxp(1), o1.Ef(1));

% perform group-BMS on data generated under the nested model
[p2, o2] = VBA_groupBMC (logEvidence_y2, options);
set (o2.options.handles.hf, 'name', 'group BMS: y_2')

fprintf('Statistics in favor of the true model (m2): pxp = %04.3f (Ef = %04.3f)\n', o2.pxp(2), o2.Ef(2));

% classical hypothesis testing
% =========================================================================
% for the sake of the example, perform same analysis using an F-test

% check if X1 better than X2 on y1
c = [zeros(2); eye(2)];
[pv,stat,df] = GLM_contrast (X1, y1, c, 'F', true);
set(gcf,'name','classical analysis of y1')

fprintf('Statistics in favor of the full model (m1): p = %04.3f (F = %04.3f)\n', pv, stat, df(1), df(2));

% check if X1 better than X2 on y2
[pv,stat,df] = GLM_contrast (X1, y2, c, 'F', true);
set(gcf,'name','classical analysis of y2')

fprintf('Statistics in favor of the full model (m1): p = %04.3f (F = %04.3f)\n', pv, stat, df(1), df(2));
% note that an absence of significance does not mean significant absence!

end

%% ########################################################################
% display subfunctions
% #########################################################################
function plotBayesFactor (logEvidence_y1, logEvidence_y2)
    [n1, x1] = VBA_empiricalDensity ((logEvidence_y1(1,:) - logEvidence_y1(2, :))');
    [n2, x2] = VBA_empiricalDensity ((logEvidence_y2(1,:) - logEvidence_y2(2 ,:))');
    hf = figure ('color' ,'w', 'name', 'demo_modelComparison: distribution of log Bayes factors');
    ha = axes ('parent', hf,'nextplot','add');
    plot (ha, x1, n1, 'color', 'r');
    plot (ha, x2, n2, 'color', 'b');
    legend (ha, {'true = model 1', 'true = model 2'});
    xlabel (ha, 'log p(y|m1) - log(y|m2)');
    ylabel (ha, 'proportion of simulations');
end
