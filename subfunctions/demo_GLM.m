% demo for GLM inference

n = 8;
X = [vec(1:n),ones(n,1)];
b = randn(2,1);
sigma = 1;
y0 = X*b;
y = y0 + randn(size(y0))/sigma;

options.inG.X = X;
g_fname = @g_GLM;
% priors.muPhi = zeros(4,1);         % prior mean on observation params
priors.SigmaPhi = 1e0*eye(2); % prior covariance on observation params
% priors.a_sigma = 1;             % Jeffrey's prior
% priors.b_sigma = 1;             % Jeffrey's prior
options.priors = priors;        % include priors in options structure
% options.inG = inG;              % input structure (grid)
dim.n_phi = 2;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states
% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);


%---- Display results ----%
displayResults(posterior,out,y,[],[],[],b,[],sigma)
