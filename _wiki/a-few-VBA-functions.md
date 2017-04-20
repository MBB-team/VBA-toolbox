---
title: "A few hidden but useful VBA stand-alone functions"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

This section lists a few VBA stand-alone functions that were designed to solve specific issues arising in the context of model-based statistical data analysis. This list is not exhaustive, and we encourage VBA users to look for (and/or contribute) other related functions within VBA's core routines.

## 1: Want to derive a sampling approximation to the states' trajectory density?

Check the function `get_MCMC_predictiveDensity.m` (and its variant `get_MCMC_predictiveDensity_fb.m` when states evolve according to some state-dependant feedback). The outputs of this procedure are time-dependant histograms, which can be eyeballed using the function `plotDensity.m`.
Note: the function `VBA_getLaplace.m` does the same thing, but under a parametric (Laplace) approximation.


## 2: Want to extract a subplot into a single matlab figure?

VBA makes intensive use of "subplots", i.e. multiple graphis within the same matlab window. These may, at times, be difficult to eyeball. Running the function `getSubplots.m` effectively attaches a "context menu" to every axis of the current matlab session. In turn, right-clicking on any axis enables to extract the subplot into a single matlab figure. This can then be used to copy-paste the graphics, e.g., into powerpoint. Note: closing the figure will automatically replace the axis at its original location.


## 3: Want to check whether any variable is numerically weird?

Numerical variables can be "weird", e.g. be complex-valued, be infinite or contain NaNs. The function `ìsweird.m` will return 1 if this is the case. Note: this function works on N-D arrays, structures, cell arrays or any combination of these.


## 4: Want to sample from an arbitrary 1D probability distribution?

Just have a look at `sampleFromarbitraryP.m` :)

> Alternatively, the function `VBA_sample.m` can be used to sample from gaussian, gamma, dirichlet, or multinomial densities...


## 5: Want to estimate impulse responses to a sequence of inputs?

Let $$x(t)$$ be the response of a system to a sequence of inputs $$u(t)$$, which can be described as a convolution operation, i.e.:

$$ x(t) = x_0 + \sum_{\tau} h\left(\tau\right) u\left(t-\tau\right) + \epsilon$$

where $$h\left(\tau\right)$$ is an unknown [impulse response function](https://en.wikipedia.org/wiki/Impulse_response) and $$\epsilon$$ are model residuals.

The function `VBA_conv2glm.m` allows you to transform a sequence of inputs into a design matrix that, when fitted using a [GLM](https://en.wikipedia.org/wiki/General_linear_model), provides an estimate of the finite impulse response function...

> The function `VBA_VolterraKernels.m` can be used to estimate the [Volterra kernels]({{ site.baseurl }}/wiki/VBA-graphical-output/) of both observed data and systems' states, given the `posterior` and `out` output VBA structures. This function exemplifies the use of `VBA_conv2glm.m`...


## 6: Want to get the square-root of a matrix?

Just have a look at `VBA_getISqrtMat.m` :)


## 7:  Want to get errorbars on estimated model residuals?

Let $$m$$ be a simple generative model of the form: $$y=g(\phi)+\epsilon$$, where $$y$$ are observed data, $$g$$ is the observation function, $$\phi$$ ar eunknown model parameters and $$\epsilon$$ are model residuals. When inverting the model (i.e. deriving the posterior density $$p\left(\theta\mid y\right)$$), VBA provides point estimates of model residuals, but does not provide the full posterior density $$p\left(\epsilon\mid y\right)$$. This can be retrieved using the function `VBA_getNoise.m`.


## 8: Deriving all n-draws from a k-urn (with replacement)

Just have a look at `VBA_getNtuples.m` :)
Note: this is problem of [combinatorics](https://en.wikipedia.org/wiki/Combinatorics), which typically arises in the context of, e.g., between-condition model comparison. This is because here, one wants to evaluate the model evidences of all possible pairings of condition-specific models...


## 9: Want to extract the states' time-dependent variances from VBA's posterior structure?

Just have a look at `VBA_getVar.m` :)


## 10: Deriving the Kullback-Leibler divergence for Normal or Gamma densities?

Just have a look at `VBA_KL.m` :)


## 11: Want to evaluate the log-evidence of the null model ($$H_0$$)?

Just have a look at `VBA_LMEH0.m` :)


## 12: Want to orthogonalize a matrix?

Sometimes, e.g. in the context of a [GLM](https://en.wikipedia.org/wiki/General_linear_model), one may want to [orthogonalize](https://en.wikipedia.org/wiki/Orthogonalization) the (design) matrix (prior to model fitting). 
This can be done using the function `VBA_orth.m`.


## 13: Want to derive the posterior probability $$P\left(\phi=0\mid y\right)$$?

Just have a look at `VBA_PP0.m` :)
Note: this function works for both evolution and observation parameters...

> The function `VBA_PPM.m` can be used to derive $$P\left(\phi\geq\mid y\right)$$ or $$P\left(a\leq\phi\leq b\mid y\right)$$. It also can be used to eyeball the corresponding posterior density...

## 14: Deriving a simple cross-validation measure of model generalizability

The predicted residual error sum of squares ([PRESS](https://en.wikipedia.org/wiki/PRESS_statistic)) statistic is a form of cross-validation used in regression analysis to provide a summary measure of the fit of a model to a sample of observations that were not themselves used to estimate the model. It is calculated as the sums of squares of the prediction residuals for those observations.
Have a look at `VBA_PRESS.m`!






