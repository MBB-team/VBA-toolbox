---
title: "VBA: results and diagnostics"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

As can be seen on [this page]({{ site.baseurl }}/wiki/VBA-graphical-output), VBA provides users with many outputs, including inversion results (e.g. parameters' estimates) and diagnostics (e.g., residuals auto-correlation function).

Below, we describe in more details VBA's main output arguments.

# Posterior density on model (unknown) variables

First of all, when calling the main inversion function, i.e..:

```matlab
[posterior, out] = ...
   VBA_NLStateSpaceModel(y, u, f_fname, g_fname, dim, options) ;
```

VBA approximates the posterior density over model variables (i.e.: evolution/observation parameters, hidden states, initial conditions and precision hyperparameters). The moments of the approximate marginal posterior densities are stored in the `posterior` structure, in a similar fashion to the `options.priors` input structure (see [this description]({{ site.baseurl }}/wiki/VBA-model-inversion-in-4-steps)), i.e.:

- **Observation parameters**
  - `posterior.muPhi`: posterior mean (vector) of observation parameters
  - `posterior.SigmaPhi`: posterior covariance (matrix) of observation parameters
- **Evolution parameters** (only for dynamical systems)
  - `posterior.muTheta`: posterior mean (vector) of evolution parameters
  - `posterior.SigmaTheta`: posterior covariance (matrix) of evolution parameters
- **Hidden states** (only for dynamical systems)
  - `posterior.muX`: posterior mean (matrix) of hidden states' time series
  - `posterior.SigmaX.current`: posterior instantaneous covariances (cell array of matrices) of hidden states
  - `posterior.SigmaX.inter`: posterior lagged covariances (cell array of matrices) of hidden states
- **Initial conditions** (only for dynamical systems)
  - `posterior.muX0`: posterior mean (vector) of initial hidden states 
  - `posterior.SigmaX0`: posterior covariance (matrix) of initial hidden states 
- **Measurement noise precision** (only for continuous data)
  - `posterior.a_sigma`: posterior shape parameter for the (Gamma-distributed) precision of observation noise
  - `posterior.b_sigma`: posterior rate parameter for the (Gamma-distributed) precision of observation noise
- **State noise precision**  (only for dynamical systems)
  - `posterior.a_alpha`: posterior shape parameter for the (Gamma-distributed) precision of state noise
  - `posterior.b_alpha`: posterior rate parameter for the (Gamma-distributed) precision of state noise

All these can be eyeballed under the 'VB inversion' tab (see [this page](VBA-graphical-output.html)).


> Posterior estimates of precision hyperparameters can be obtained from their posterior mean, which reduces to the ratio of shape and rate parameters, e.g. (for the observation noise $$\sigma$$):
$$E[\sigma|y]= \frac{a_{\sigma}}{b_{\sigma}}$$
where $$a_{\sigma}$$ and $$b_{\sigma}$$ are the posterior shape and rate parameters, respectively.


# Model quality metrics

The above inversion call also returns the structure `out`, which can be queried for all sorts of diagnostics. Among these, VBA computes model quality metrics:

- **Model's log-evidence**: `out.F`. This is the variational Bayesian approximation to the model's marginal likelihood, which is required for model comparison. A dummy "Bayesian p-value" for the model can then be computed as follows:
   
  ```matlab
dF = out.F - out.diagnostics.LLH0 ;
bayesianP = 1./(1+exp(dF)) ;
```

where `out.diagnostics.LLH0` is the evidence for the "null" (i.e. a model that assumes that data are random samples). The above "Bayesian p-value" is in fact the posterior probability of the null...
  
> More generally, the model's log-evidence is used for model comparison purposes. This specific issue is treated in other sections. For example, group-level model selection is described [here]({{ site.baseurl }}/wiki/BMS-for-group-studies), whereas model selection with large model spaces is described [here]({{ site.baseurl }}/wiki/Comparing-large-spaces-of-models).
  
- **Goodness-of-fit metrics**:
  - percentage of variance explained: `out.fit.R2`.
  - log-likelihood: `out.fit.LL`.
  - [AIC](https://en.wikipedia.org/wiki/Akaike_information_criterion): `out.fit.AIC`
  - [BIC](https://en.wikipedia.org/wiki/Bayesian_information_criterion): `out.fit.BIC`
  - when dealing with categorical (binary) data, VBA also outputs the classification accuracy, i.e. the percentage of correct classifications (see `out.fit.acc`), as well as the so-called ["balanced classification accuracy"](http://www.sciencedirect.com/science/article/pii/S1053811913002371) (see `out.fit.bacc`). 

> **TIP**: here AIC and BIC scores are defined as the log-likelihood, minus a model complexity penalty term (typically related to the number of unknown model variables). This means that for all model quality metrics, the higher the score, the better the model.
All these are given under the `summary` tab (see [this page]({{ site.baseurl }}/wiki/VBA-graphical-output)). Note that all goodness-of-fit metrics can be re-derived using the function `VBA_fit.m` (which simply requires both `posterior` and `out` structures).

# Inversion diagnostics

Typically, four types of diagnostics can be eyeballed in a systematic manner:

- **model parameters' posterior correlation matrix**: `out.diagnostics.C`. This matrix is useful for checking potential [non-identifiability](https://en.wikipedia.org/wiki/Identifiability) issues, which would express as non-zero posterior correlation between model parameters.
- **residuals empirical auto-correlation**: `out.diagnostics.dy.R`. This is useful for checking the absence of structure in model residuals, which would signal "underfitting".
- **residuals empirical histogram**: this can be used as a post-hoc sanity check for a key prior assumption, under which the likelihood function is derived, namely: residuals should be normally distributed. Note: although, formally speaking, this does not apply to binary data, a well-behaved model fit under a binomial likelihood would exhibit a normal distribution of empirical residuals.
- **micro-time resolution hidden states and predicted data**: this can be used to eyeball the model predictions interpolated outside the sampling grid. Note: this only applies to models that rely on (almost) continuous dynamical systems, whereby the evolution function is applied more than one time in between two time samples.

These diagnostics can be eyeballed under the 'diagnostics' tab (see [this page]({{ site.baseurl }}/wiki/VBA-graphical-output)).

In addition, for dynamical systems, VBA provides:

- **Volterra kernels**: `out.diagnostics.kernels`. This diagnostic analysis allows one to identify the system's states' and observables' impulse response to experimentally controlled inputs (even when the system is strongly nonlinear). Volterra kernels can be eyeballed under the 'kernels' tab (see [this page]({{ site.baseurl }}/wiki/VBA-graphical-output)).
