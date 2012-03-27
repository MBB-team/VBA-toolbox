% Demo for Savage-Dickey nested model comparison
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

% close all
clear variables

% Choose basic settings for simulations
sigma = 1e-1;            % precision 
g_fname = @g_GLM;        % observation function


% Build priors structure
priors.muPhi = zeros(3,1);         % prior mean on observation params
priors.SigmaPhi = 1e4*eye(3); % prior covariance on observation params
priors.a_sigma = 1e1;             % Jeffrey's prior
priors.b_sigma = 1e1;             % Jeffrey's prior
options.priors = priors;        % include priors in options structure
options.DisplayWin = 1;
options.verbose = 0;
options.Laplace = 1;
options.MinIter = 30;
options.updateHP = 1;
dim.n_phi = 3;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states

% fill in priors
[options,u,dim] = VBA_check(zeros(1e2,1),[],[],g_fname,dim,options);
priors = options.priors;
priors2 = priors;
priors2.SigmaPhi(2,2) = 0;

N = 1e1;

F1 = zeros(N,2);
F2 = zeros(N,2);
F2sd = zeros(N,2);

F1_0 = zeros(N,2);
F2_0 = zeros(N,2);
F2sd_0 = zeros(N,2);

for i=1:N

    i
    
    % simulate data under full and reduced models
    phi = 4+ randn(3,1);
    inG.X = [randn(1e2,2),ones(1e2,1)];
    [gx] = feval(g_fname,[],phi,[],inG);
    y1 = gx + sqrt(sigma.^-1)*randn(size(gx));
    phi(2) = 0;
    [gx] = feval(g_fname,[],phi,[],inG);
    y2 = gx + sqrt(sigma.^-1)*randn(size(gx));
    options.inG = inG; 

    % Invert full model on 'full' data
    options.priors = priors;
    [p1,o1] = VBA_NLStateSpaceModel(y1,[],[],g_fname,dim,options);
    F1(i,1) = o1.F;
    
    % Invert reduced model on 'full' data
    options.priors = priors2;
    [p2,o2] = VBA_NLStateSpaceModel(y1,[],[],g_fname,dim,options);
    F2(i,1) = o2.F;
    
    % Use Savage-Dickey ratio
    [F2sd(i,1),p2sd] = VB_SavageDickey(p1,priors,o1.F,dim,priors2);
    
    % Use frequentist limit:
    F1_f(i,1) = lev_GLM(y1,inG.X);
    F2_f(i,1) = lev_GLM(y1,inG.X(:,[1,3]));
    
    % repeat without Laplace free energy:
    [F1_0(i,1)] = getF(p1,o1,0);
    [F2_0(i,1)] = getF(p2,o2,0);
    [F2sd_0(i,1)] = VB_SavageDickey(p1,priors,F1_0(i,1),dim,priors2);
    
    
    % Invert full model on 'reduced' data
    options.priors = priors;
    [p1,o1] = VBA_NLStateSpaceModel(y2,[],[],g_fname,dim,options);
    F1(i,2) = o1.F;
    
    % Invert reduced model on 'reduced' data
    options.priors = priors2;
    [p2,o2] = VBA_NLStateSpaceModel(y2,[],[],g_fname,dim,options);
    F2(i,2) = o2.F;
    
    % Use Savage-Dickey ratio
    [F2sd(i,2),p2sd] = VB_SavageDickey(p1,priors,o1.F,dim,priors2);
    
    % Use frequentist limit:
    F1_f(i,2) = lev_GLM(y2,inG.X);
    F2_f(i,2) = lev_GLM(y2,inG.X(:,[1,3]));
    
    % repeat without Laplace free energy:
    [F1_0(i,2)] = getF(p1,o1,0);
    [F2_0(i,2)] = getF(p2,o2,0);
    [F2sd_0(i,2)] = VB_SavageDickey(p1,priors,F1_0(i,2),dim,priors2);
    
    

end

hf = figure('color',[1 1 1]);
subplot(2,2,1),plot(F1-F2sd,F1-F2,'.')
xlabel('log BF: Savage-Dickey')
ylabel('log BF: full inversion')
legend({'log p(m_f|y) - log p(m_r|y): Y~p(y|m_f)','log p(m_f|y) - log p(m_r|y): Y~p(y|m_r)'})
title('under the Laplace approx')
grid on

subplot(2,2,2),plot(F1_0-F2sd_0,F1_0-F2_0,'.')
xlabel('log BF: Savage-Dickey')
ylabel('log BF: full inversion')
legend({'log p(m_f|y) - log p(m_r|y): Y~p(y|m_f)','log p(m_f|y) - log p(m_r|y): Y~p(y|m_r)'})
title('without the Laplace approx')
grid on


subplot(2,2,3),plot((F1-F2)-(F1_0-F2_0),'.')
xlabel('simulations')
ylabel('log BF: Laplace - no Laplace')
legend({'Y~p(y|m_f)','Y~p(y|m_r)'})
title('with and wihout Laplace')
grid on


subplot(2,2,4),plot(F1_f-F2_f,'.')
xlabel('simulations')
ylabel('log BF: frequentist limit')
legend({'Y~p(y|m_f)','Y~p(y|m_r)'})
title('frequentist limit')
grid on


