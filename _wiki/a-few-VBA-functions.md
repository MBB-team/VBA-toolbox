---
title: "A few useful VBA stand-alone functions"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

This section lists a few VBA stand-alone functions that were designed to solve specific issues arising in the context of model-based statistical data analysis. This list is not exhaustive, and we encourage VBA users to look for (and/or contribute) other related functions within VBA's core routines.

## Want to derive a sampling approximation to the states' trajectory density?

Check the function `get_MCMC_predictiveDensity.m` (and its variant `get_MCMC_predictiveDensity_fb.m` when states evolve according to some state-dependant feedback). The outputs of this procedure are time-dependant histograms, which can be eyeballed using the function `plotDensity.m`.
Note: the function `VBA_getLaplace.m` does the same thing, but under a parametric (Laplace) approximation.


## Want to check whether any variable is numerically weird?

Numerical variables can be "weird", e.g. be complex-valued, be infinite or contain NaNs. The function `VBA_isweird.m` will return 1 if this is the case. Note: this function works on N-D arrays, structures, cell arrays or any combination of these.


## Want to sample from an arbitrary 1D probability distribution?

Have a look at `VBA_random('Arbirary', p, values)`

> Alternatively, the function `VBA_random.m` can be used to sample from gaussian, gamma, dirichlet, or multinomial densities...


## Want to estimate impulse responses to a sequence of inputs?

Let $$x(t)$$ be the response of a system to a sequence of inputs $$u(t)$$, which can be described as a [convolution](https://en.wikipedia.org/wiki/Convolution) operation, i.e.:

$$ x(t) = x_0 + \sum_{\tau} h\left(\tau\right) u\left(t-\tau\right) + \epsilon$$

where $$h\left(\tau\right)$$ is an unknown [impulse response function](https://en.wikipedia.org/wiki/Impulse_response), $$\tau$$ is the convolution lag and $$\epsilon$$ are model residuals.

The function `VBA_conv2glm.m` allows you to transform a sequence of inputs into a design matrix that, when fitted using a [GLM](https://en.wikipedia.org/wiki/General_linear_model), provides an estimate of the finite impulse response function...

> The function `VBA_VolterraKernels.m` can be used to estimate the [Volterra kernels]({{ site.baseurl }}/wiki/Volterra-decomposition/) of both observed data and systems' states, given the `posterior` and `out` output VBA structures. This function exemplifies the use of `VBA_conv2glm.m`...


## Want to get the square-root of a matrix?

Have a look at `VBA_sqrtm.m` :)


## Want to get errorbars on estimated model residuals?

Let $$m$$ be a simple generative model of the form: $$y=g(\phi)+\epsilon$$, where $$y$$ are observed data, $$g$$ is the observation function, $$\phi$$ are unknown model parameters and $$\epsilon$$ are model residuals. When inverting the model (i.e. deriving the posterior density $$p\left(\theta\mid y\right)$$), VBA provides point estimates of model residuals, but does not provide the full posterior density $$p\left(\epsilon\mid y\right)$$. This can be retrieved using the function `VBA_getNoise.m`.


## Deriving all n-draws from a k-urn (with replacement)

Have a look at `VBA_getNtuples.m` :)

> This is problem of [combinatorics](https://en.wikipedia.org/wiki/Combinatorics), which arises, e.g., in the context of [between-condition model comparison]({{ site.baseurl }}/wiki/BMS-for-group-studies). This is because here, one wants to evaluate the model evidences of all possible pairings (n-tuples, where $$n$$ is the number of conditions) of condition-specific models (where there are $$k$$ models)...


## Want to extract the states' variances from VBA's posterior structure?

Have a look at `VBA_getVar.m` :)


## Deriving the Kullback-Leibler divergence for Normal or Gamma densities

Have a look at `VBA_KL.m` :)


## Want to evaluate the log-evidence of the null model?

Have a look at `VBA_LMEH0.m` :)

> The 'null model' here is such that the data is entirely random. It may or may not directly correspond to a given 'null assumption' ($$H_0$$) used in classical testing procedures.


## Want to orthogonalize a matrix?

Have a look at `VBA_orth.m` :)


## Want to derive the posterior probability $$P\left(\phi=0\mid y\right)$$?

Have a look at `VBA_PP0.m` :)

Note: this function works for both evolution and observation parameters...

> The function `VBA_PPM.m` can be used to derive $$P\left(\phi\geq a\mid y\right)$$ or $$P\left(a\leq\phi\leq b\mid y\right)$$. It also can be used to eyeball the corresponding posterior density...


## Deriving a simple cross-validation measure of model generalizability

Have a look at `VBA_PRESS.m` :)

> The predicted residual error sum of squares ([PRESS](https://en.wikipedia.org/wiki/PRESS_statistic)) statistic is a form of cross-validation used in regression analysis to provide a summary measure of the fit of a model to a sample of observations that were not themselves used to estimate the model. It is calculated as the sums of squares of the prediction residuals for those observations.


## Want to derive Savage-Dickey ratios?

Have a look at `VBA_SavageDickey.m` :)

> [Savage-Dickey ratios](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0059655) can be used to compute the free energy and posterior moments of a reduced model, given the prior and posterior densities of a full model. The key added-value of this scheme is that one does not need to perform the corresponding model inversion! This is particularly useful when [comparing models within large model spaces]({{ site.baseurl }}/wiki/Comparing-large-spaces-of-models).


## Tools for classical contrast testing.

VBA includes a function (`GLM_contrast.m`) that enables [classical (i.e. frequentist) statistical testing]({{ site.baseurl }}/wiki/statistical-models) in the context of a GLM. Below is a set of related useful functions:

- `Contrast_MEbins.m`: derives the [contrast matrix](https://en.wikipedia.org/wiki/Contrast_(statistics)) for an F-test of group mean differences (useful when the number of groups or conditions is bigger than 2).
- `FtoR2.m`: converts an F-statistic into a [coefficient of determination](https://en.wikipedia.org/wiki/Coefficient_of_determination) (R2).
- `findCI.m`: finds [confidence intervals](https://en.wikipedia.org/wiki/Confidence_interval) for given t or F statistics.
- `GLM_tolerance.m`: computes [regression tolerance](https://en.wikipedia.org/wiki/Multicollinearity) of a given design matrix.
- `lev_GLM.m`: computes the log-evidence of a GLM (frequentist limit).
- `PRESS_GLM.m`: evaluates the [PRESS-R2](https://en.wikipedia.org/wiki/PRESS_statistic) cross-validation metric for a GLM.
- `removeOutliers.m`: detects [outliers](https://en.wikipedia.org/wiki/Outlier) based upon robust moment-matched Gaussian distribution.
- `testPower.m`: returns the [statistical power](https://en.wikipedia.org/wiki/Statistical_power) of a test, given the expected [effect size](https://en.wikipedia.org/wiki/Effect_size) and the corresponding degrees-of-freedom.
- `GLM_2pieces.m`: fits a [piece-wise linear model](https://en.wikipedia.org/wiki/Piecewise_linear_function).


## Want to extract a subplot into a single matlab figure?

VBA makes intensive use of "subplots", i.e. multiple graphis within the same matlab window. These may, at times, be difficult to eyeball. Running the function `VBA_getSubplots.m` effectively attaches a "context menu" to every axis of the current matlab session. In turn, right-clicking on any axis enables one to extract the subplot into a single matlab figure. This can then be used to copy-paste the graphics, e.g., into powerpoint. Note: closing the figure will automatically replace the axis at its original location.
