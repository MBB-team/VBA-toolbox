---
title: "Results and diagnostics"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

As can be seen on [this page]({{ site.baseurl }}/wiki/VBA-graphical-output), VBA provides users with many outputs, including inversion results (e.g. parameters' estimates) and diagnostics (e.g., residuals auto-correlation function).

Below, we describe in more details VBA's main output arguments.

# Posterior density on model (unknown) variables

First of all, VBA derives an approximation to the posterior density over model variables (i.e.: evolution/observation parameters, hidden states, initial conditions and precision hyperparameters). When calling the main inversion function, e.g.:

```matlab
[posterior, out] = ...
   VBA_NLStateSpaceModel(y, u, f_fname, g_fname, dim, options) ;
```

moments of the approximate posterior density are stored in `posterior`, in a similar fashion than the `options.priors` input structure (see [this description]({{ site.baseurl }}/wiki/VBA-model-inversion-in-4-steps)), i.e.:

- **Observation parameters**
  - `posterior.muPhi`
  - `posterior.SigmaPhi`
- **Evolution parameters** (only for dynamical systems)
  - `posterior.muTheta`
  - `posterior.SigmaTheta`
- **Hidden states** (only for dynamical systems)
  - `posterior.muX`
  - `posterior.SigmaX.current` (a cell array of instantaneous states' covariance matrices)
  - `posterior.SigmaX.inter` (a cell array of lagged- covariance matrices)
- **Initial conditions** (only for dynamical systems)
  - `posterior.muX0`
  - `posterior.SigmaX0`
- **Measurement noise precision** (only for continuous data)
  - `posterior.a_sigma`
  - `posterior.b_sigma`
- **State noise precision**  (only for dynamical systems)
  - `posterior.a_alpha`
  - `posterior.b_alpha`

All these can be eyeballed under the 'VB inversion' tab (see [this page](VBA-graphical-output.html)).

# Model quality metrics

The above inversion call also returns the structure `out`, which can be queried for all sorts of diagnostics. Among these, VBA computes model quality metrics:

- **Model's log-evidence**: `out.F`. This is the variational Bayesian approximation to the model's marginal likelihood, which is required for model comparison. A dummy "Bayesian p-value" for the model can then be computed as follows:
   
  ```matlab
dF = out.F - out.diagnostics.LLH0 ;
bayesianP = 1./(1+exp(dF)) ;
```

  where `out.diagnostics.LLH0` is the evidence for the "null" (i.e. a model that assumes that data are random samples). The above "Bayesian p-value" is in fact the posterior probability of the null...
- **Fit accuracy metrics**:
  - percentage of variance explained: `out.fit.R2`. NB: when dealing with categorical (binary) data, `out.fit.R2` is the balanced classification accuracy, i.e. the percentage of correct classifications (based upon a 0.5 threshold).
  - log-likelihood: `out.fit.LL`.
  - AIC: `out.fit.AIC`
  - BIC: `out.fit.BIC`

> **TIP**: here AIC and BIC scores are defined as the log-likelihood, minus a model complexity penalty term (typically related to the number of unknown model variables). This means that for all model quality metrics, the higher the score, the better the model.
All these are given under the 'summary' tab (see [this page]({{ site.baseurl }}/wiki/VBA-graphical-output)).

# Inversion diagnostics

Typically, three types of diagnostics can be eyeballed in a systematic manner:

- **model parameters' posterior correlation matrix**: `out.diagnostics.C`. This matrix is useful for checking potential non-identifiability issues, which would express as non-zero posterior correlation between model parameters.
- **residuals empirical auto-correlation**: `out.diagnostics.dy.R`. This is useful for checking the absence of structure in model residuals, which would signal "underfitting".
- **residuals empirical histogram**: this can be used as a post-hoc sanity check for a key prior assumption, under which the likelihood function is derived, namely: residuals should be normally distributed. Note: although, formally speaking, this does not apply to binary data, a well-behaved model fit under a binomial likelihood would exhibit a normal distribution of empirical residuals.
- **micro-time resolution hidden states and predicted data**: this can be used to eyeball the model predictions interpolated outside the sampling grid. Note: this only applies to models that rely on (almost) continuous dynamical systems, whereby the evolution function is applied more than one time in between two time samples.

These diagnostics can be eyeballed under the 'diagnostics' tab (see [this page]({{ site.baseurl }}/wiki/VBA-graphical-output)).

- **Volterra kernels**: `out.diagnostics.kernels`. This diagnostic analysis allows one to identify the system's states' and observables' impulse response to experimentally controlled inputs (even when the system is strongly nonlinear). Volterra kernels can be eyeballed under the 'kernels' tab (see [this page]({{ site.baseurl }}/wiki/VBA-graphical-output)).
