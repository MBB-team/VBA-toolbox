% demo for binomial data inversion
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

clear variables
% close all


% p = 2e2;
phi = [-1;0]; % phi(1) = log sigmoid slope , phi(2) = inflexion point
% u = 1e0*randn(p,1);
gridu = -10:1e-1:10;
u = gridu(:);
p = length(u);
% u = sort(u);

% get binomial sampling probabilities
sx = g_sigm_binomial([],phi,u,[]);

% sample binomial data
y = zeros(p,1);
seed = 1e4*rand;
for t=1:p
    try
        [y(t),seed] = binomial_sample(1,sx(t),seed);
    catch
        [y(t)] = sampleFromArbitraryP([sx(t),1-sx(t)],[1,0],1);
    end
end

figure
plot(sx,'r')
hold on
plot(y,'k.')
legend({'sampling proba: p(y|x)','data samples: y'})

dim.n_phi = 2;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states
dim.n_t = 1;
dim.p = p;

g_fname = @g_sigm_binomial;     % observation function

% Call inversion routine
options.binomial = 1;
options.gradF = 0;
options.priors.SigmaPhi = 1e-1*eye(2);
options.DisplayWin = 0;
options.verbose = 0;

try
    [e] = designEfficiency([],g_fname,dim,options,u,'parameters');
    disp(['Design efficiency: ',num2str(e)])
catch
    disp('Warning: could not estimate design efficiency!')
end


mu = zeros(dim.n_phi,p);
va = zeros(dim.n_phi,p);
hf = figure;
ha = gca(hf);
set(ha,'nextplot','add')
for t=1:p
    dim.p = t;
    [posterior,out] = VBA_NLStateSpaceModel(...
        y(1:t),u(1:t),[],g_fname,dim,options);
    mu(:,t) = posterior.muPhi;
    va(:,t) = diag(posterior.SigmaPhi);
    if t > 1
        cla(ha)
        [haf,hf] = plotUncertainTimeSeries(mu(:,1:t),sqrt(va(:,1:t)),1:t,ha,1:2);
    end
end


%---- Display results ----%
displayResults(posterior,out,y,[],[],[],phi,[],[])
