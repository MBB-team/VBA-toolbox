function results = VBA_sensitivityAnalysis(posterior,out,filter,N_RUNS)
%% 
% VBA_SENSIVITYANALYSIS Given a behavioural DCM, compute the sensitivity of
% the behavioural responses to the functionnal connectivity parameters.
%
% results = VBA_sensitivityAnalysis(posterior,out,[N_RUNS=64*nConnections])
%
%

warning off


% select only connexion of interest
% -------------------------------------------------------------------------

inF= out.options.inF;
Aself = inF.A;
Aself(Aself==1) = inF.indA;
true_As = setdiff(inF.indA,diag(Aself));
idxTheta=[true_As inF.indB{:} inF.indC inF.indD{:}]; %inF.indC


if ~isempty(filter)
    idxTheta = intersect(idxTheta,find(diag(posterior.SigmaTheta)>1e-6));
    idxTheta = intersect(idxTheta,find(filter));
end

thetas_lbl={};
for i=1:length(true_As)
    thetas_lbl{end+1} = sprintf('A^ _%d',i);
end
for k=1:length(inF.indB)
    for i=1:length(inF.indB{k})
    thetas_lbl{end+1} = sprintf('B^%d_%d',k,i);
    end
end
for i=1:length(inF.indC)
    thetas_lbl{end+1} = sprintf('C_%d',i);
end
for k=1:length(inF.indD)
    for i=1:length(inF.indD{k})
    thetas_lbl{end+1} = sprintf('D^%d_%d',k,i);
    end
end

%% compute the minimum number of runs
if nargin<4
    N_RUNS = numel(idxTheta) * 3 * 50 ;
end
fprintf('running %d simulations:\n',N_RUNS);


%% Run simulations over connectivity space and compute corresponding kernels
resps = [out.options.sources(2:end).out];
nResps = numel(resps);
n_t = out.dim.n_t;
n_u = out.dim.u;

% runArray = zeros(N_RUNS,nResps,n_t);
thetaArray = zeros(N_RUNS,numel(idxTheta));
kernel(N_RUNS,nResps) = struct('timeline',[],'params',[],'timeseries',[],'landmarks',struct(),'sigma',struct());


kout=[];
parfor iRun=1:N_RUNS
    % get a new set of connections
    [posterior_run,thetaArray(iRun,:)] = perturb_theta(posterior,idxTheta);
    % simulate a run and compute the volterra kernels
    try 
        kernel(iRun,:) = find_kernel(posterior_run,out); % ,runArray(iRun,:,:)
    catch e
        e.message
        kout = [kout iRun];
    end
    fprintf('.');
end
kout;
kernel(kout,:) = [];
% runArray(kout,:,:) = [];
thetaArray(kout,:) = [];

%% compute regressors
thetaArray_z = zscore(thetaArray) ;

% compute relative weights of each parameter on each time step/observation
% -------------------------------------------------------------------------
% betaArray = zeros(numel(idxTheta)+1,nResps,n_t);
% for iObs = 1:nResps
%     parfor t=1:n_t  
%         betaArray(:,iObs,t) = glmfit(thetaArray_z,runArray(:,iObs,t),'normal');
%     end
% end

% compute effect of parameters on kernels landmarks
% -------------------------------------------------------------------------
kernelLandmarksBeta = zeros(nResps,n_u,1+numel(idxTheta));
kernelLandmarksBetaSe = zeros(nResps,n_u,1+numel(idxTheta));
for iObs = 1:nResps
    kernel_obs = kernel(:,iObs);
    for iu = 1:n_u   
        aMax = arrayfun(@(k) k.landmarks(iu).aMax,kernel_obs);
        aMaxSe = arrayfun(@(k) k.sigma.landmarks(iu),kernel_obs);
        [beta,stats]=robustfit(thetaArray_z,aMax);
        kernelLandmarksBeta(iObs,iu,:) = beta;
        kernelLandmarksBetaSe(iObs,iu,:) = stats.s *sqrt(N_RUNS) ;
    end
end


%%
% Store results
% -------------------------------------------------------------------------
results.posterior=posterior;
results.out=out;
results.thetaArray = thetaArray;
results.theta.lbl = thetas_lbl;
results.theta.idx = idxTheta;
% results.timeseriesBeta = betaArray(2:end,:,:);
results.kernel = kernel;
results.kernelLandmarksBeta = kernelLandmarksBeta;
results.kernelLandmarksBetaSe = kernelLandmarksBetaSe;

%%
% Display results
% -------------------------------------------------------------------------
try
    VBA_sensitivityAnalysisDisplay(results)
end
end

function [posterior,perturbed_theta]=perturb_theta(posterior,idx)
    sdTheta = sqrt(diag(posterior.SigmaTheta(idx,idx)));
    perturbed_theta = posterior.muTheta(idx) + sdTheta.*randn(numel(idx),1);
    posterior.muTheta(idx) = perturbed_theta;
end

% function [posterior,perturbed_phi]=perturb_phi(posterior,idx)
%     varPhi = sqrt(diag(posterior.SigmaPhi(idx,idx)));
%     perturbed_phi = posterior.muPhi(idx) + varPhi.*randn(numel(idx),1);
%     posterior.muPhi(idx) = perturbed_phi;
% end

function beta = getBeta(thetas, values)

for iVal=1:size(values,2);
    ns = 2;
    g_fname = @g_GLM;

    % Build options structure
    inG.X = [ones(size(theta,1),1) thetas];
    options.inG     = inG;
    priors.a_sigma = 1e0;
    priors.b_sigma = 1e0;
    priors.muPhi = zeros(ns,1);
    priors.SigmaPhi = 1e4*eye(ns);
    options.priors      = priors;
    dim.n_theta         = 0;
    dim.n_phi           = ns;
    dim.n               = 0;


    % Invert model with confounds
    [p,~] = VBA_NLStateSpaceModel(values(:,iVal),[],[],g_fname,dim,options);
    
    beta.mu(:,iVal) = p.muPhi;
    beta.sigma(:,iVal) = diag(p.SigmaPhi);
end

end

