---
title: "A few hidden but useful VBA stand-alone functions"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

This section lists a few VBA stand-alone functions that were designed to solve specific issues arising in the context of model-based statistical data analysis. This list is not exhaustive, and we encourage VBA users to look for (and/or contribute) other related functions within VBA's core routines.

# 1: Want to derive a sampling approximation to the states' trajectory density?

Check the function `get_MCMC_predictiveDensity.m` (and its variant `get_MCMC_predictiveDensity_fb.m` when states evolve according to some state-dependant feedback). The outputs of this procedure are time-dependant histograms, which can be eyeballed using the function `plotDensity.m`.
Note: the function `VBA_getLaplace.m` does the same thing, but under a parametric (Laplace) approximation.


# 2: Want to extract a subplot into a single matlab figure?

VBA makes intensive use of "subplots", i.e. multiple graphis within the same matlab window. These may, at times, be difficult to eyeball. Running the function `getSubplots.m` effectively attaches a "context menu" to every axis of the current matlab session. In turn, right-clicking on any axis enables to extract the subplot into a single matlab figure. This can then be used to copy-paste the graphics, e.g., into powerpoint. Note: closing the figure will automatically replace the axis at its original location.

# 3: Want to check whether any variable is numerically weird?

Numerical variables can be "weird", e.g. be complex-valued, be infinite or contain NaNs. The function `ìsweird.m` will return 1 if this is the case. Note: this function works on N-D arrays, structures, cell arrays or any combination of these.


# 4: Want to sample from an arbitrary 1D probability distribution?

Just have a look at `sampleFromarbitraryP.m` :)

# 5: Want to estimate finite impulse responses to a sequence of inputs?

Let $$x(t)$$ be the response of a system to a sequence of inputs $$u(t)$$, which can be described as a convolution operation, i.e.:

$$ x(t) = x_0 + \sum_{\tau} h\left(\tau\right) u\left(t-tau)\right) + \epsilon$$

where $$h\left(\tau\right)$$ is an unknown finite [impulse response function](https://en.wikipedia.org/wiki/Impulse_response).







