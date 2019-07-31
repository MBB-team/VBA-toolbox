---
title: "Model identifiability and confusion analyses"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


From a statistical perspective, model-based data analysis reduces to either estimating unknown model parameters or comparing candidate models (given experimental data). However, such analyses may not always be accurate. Below we highlight numerical analyses that can be performed to address the question of whether parameter estimation and/or model selection are indeed reliable, under a given set of experimental constraints (cf. design, data quality, etc...). 

 
## Simulation-recovery analysis

Model *identifiability* analysis aims at answering a central question: can parameters be identified from observed data? The answer to this question is not trivial. First of all, it typically depends upon the experimental design. This implies that, ideally, one should make sure the design is compatible with the ensuing model-based data analysis. Second, it depends upon the signal-to-noise ratio. This is because noise may either mask informative variations in the data, or be confused with data features that would otherwise be caused by the model. Third, it depends upon the generative model. For example, two parameters may have a similar impact on the data. This would cause some identifiability issue...

This is why performing an identifiability analysis is always a healthy counterpart to any model-based data analysis. In what follows, we sketch what we call a **simulation-recovery** analysis, which is one out of many ways for assessing model identifiability:

```
1) for i=1:N (Monte-Carlo simulations)
      sample model parameters under the prior distribution
      simulate data given simulated parameters
      invert model on simulated data and store estimated parameters
   end
2) regress the estimated parameters on simulated parameters
```

Of particular interest here is the relative amount of variance in each estimated parameter that can be explained by variations in simulated parameters. In fact, any non-diagonal element in the matrix of regression coefficients signals a potential non-identifiability issue between the corresponding parameters. Note that the total amount of variance explained by the simulated parameters is also of interest. This is because a low amount of explained variance means that the identifiability of the corresponding parameter depends upon specific combinations of other parameters (high-order interactions). This may arise in the context of strong nonlinearities in the generative model...  


> TIP: The simulated parameters can be sampled under the prior distributions that are used for model inversion. Alternatively, one can sample them under their empirical distribution (this can be done after the data analysis has been performed).   


## Model confusion analysis

Recall that a given generative model is specified in terms of observation/evolution functions, as well as priors on model parameters. Obviously, potential causes of model non-identifiabilty (non-informative design, low signal-to-noise ratio, parameter redundancy, etc...) are also issues for model selection. In what follows, we describe a simple variant of **confusion** analysis, which aims at assessing whether distinct models may be confused with each other (given the available experimental data):

```
1) for i=1:N (Monte-Carlo simulations)
      for sm=1:M [loop over simulated models]
          simulate data under model "sm"
          for cm=1:M [loop over candidate models]
              invert model "cm" on simulated data
          end
          perform bayesian model selection
       end
   end
2) confusion matrix = frequency with which each candidate model is selected (for each simulated model)
```

Here again, any non-diagonal element in the **confusion matrix** signals a potential confusion between the selected model and the true (hidden) model...

 > TIP: The accuracy of such Monte-Carlo simulations depends upon the number `N` of simulations performed. If time permits, we suggest to set `N` at least at two orders of magnitude greater than the number of parameters (identifiability analyses) or models (confusion analyses).
