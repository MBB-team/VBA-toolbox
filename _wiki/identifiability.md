---
title: "Model identifiability and confusion analyses"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


 
## Simulation-recovery analysis

Model *identifiability* analysis aims at answering a central question: can parameters be identified from observed data? The answer to this question is not trivial. First of all, it typically depends upon the experimental design. This implies that, ideally, one should make sure the design is compatible with the ensuing model-based data analysis. Second, it depends upon the signal-to-noise ratio. This is because noise may either mask informative variations in the data, or be confused with data features that would otherwise be caused by the model. Third, it depends upon the generative model. For example, two parameters may have a similar impact on the data. This would cause some identifiability issue...

This is why performing an identifiability analysis is alsways a healthy counterpart to any model-based data analysis. In what follows, we describe what we call a **simulation-recovery** analysis, which is one out of many ways for assessing model identifiability:

```
1) for i=1:N (Monte-Carlo simulations)
      sample model parameters under the prior distribution
      simulate data given simulated parameters
      invert model on simulated data and store estimated parameters
   end
2) regress the estimated parameters on simulated parameters
```

Of particular interest here is the relative amount of variance in each estimated parameter that can be explained by variations in simulated parameters.

Below, we provide concrete example on a modified Q-learning model, which has 2 observation parameters (inverse temperature and bias), 2 evolution parameters (learning rate and sensitivity to reward) and 2 hidden states (the values of each option).

```
% set generative model
nt = 50; % number of trials in the learning task
f_fname = @f_Qlearn3; % evolution function (Q-learning)
g_fname = @g_softmax; % observation function (softmax mapping)
options.binomial = 1;
options.skipf = zeros(1,nt);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

% set feedback structure (cf. Q-learning model)
fb.inH.frame = 'gain'; % could be 'loss'
fb.inH.Rsize = 1; % feedback magnitude
fb.inH.f = [0.8;0.2]; % action-specific gain frequency
fb.h_fname = @h_gainsAndLosses;
fb.indy = 1;
fb.indfb = 2;

% Monte-Carlo simulations
Nmc = 5e2; % number of Monte-Carlo simulations
Y = NaN(Nmc,6); % this will be used to store estimated params
X = X(Nmc,6);  % this will be used to store simulated params
for imc = 1:Nmc
    % simulate data
    theta = randn(2,1);
    phi = randn(2,1);
    x0 = rand(2,1);
    [y,x,x0,eta,e,u] = simulateNLSS_fb(nt,f_fname,g_fname,theta,phi,zeros(2,nt),Inf,Inf,options,x0,fb);
    % invert model on simulated data
    options.DisplayWin = 0;
    dim = struct('n',2,'n_theta',2,'n_phi',2);
    [posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
    % store simulated and estimated parameters
    Y(imc,:) = [posterior.muPhi',posterior.muTheta',posterior.muX0'];
    X(imc,:) = [phi',theta',x0'];
end

% regress estimated params on simulated params
zY = zscore(Y);
zX = zscore(X);
[pv,stat,df,all] = GLM_contrast(zX,zY,eye(6),'F',1);
b = all.b; % regression coefficients
R = all.R2; % amount of explained variance
figure,imagesc(b)
```

Any non-diagonal element in the matrix of regression coefficients signals a potential non-identifiability issue between the corresponding parameters. Note that the total amount of variance explained by the simulated parameters is also of interest. This is because a low amount of explained variance means that the identifiability of the corresponding parameter depends upon specific combinations of other parameters (high-order interactions). This may arise in the context of strong nonlinearities in the generative model...  


> TIP: Here, the simulated parameters were sampled under the prior distributions that were used for model inversion. Alternatively, one can sample them under their empirical distribution (this can be done after the data analysis has been performed).   


## Model confusion analysis

Recall that a given generative model is specified in terms of observation/evolution functions, as well as priors on model parameters. Obviously, potential causes of model non-identifiabilty (non-informative design, low signal-to-noise ratio, parameter redundancy, etc...) are also issues for model selection. In what follows, we describe a simple variant of **confusion** analysis, which aims at assessing whether distinct models may be confused with each other (given the available experimental data):

```
1) for i=1:N (Monte-Carlo simulations)
      for sm=1:M [loop over simulated models]
          simulate data under simulated model 'sm'
          for im=1:M [loop over estimated models]
              invert model 'im' on simulated data
          end
          perform bayesian model selection
       end
   end
2) confusion matrix = frequency with which each candidate model is selected (for each simulated model)
```

Here again, any non-diagonal element in the **confusion matrix** signals a potential confusion between the selected model and the true (hidden) model...

 
