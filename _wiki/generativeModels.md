---
title: "How general is the class of VBA's generative models?"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

As we have described [here]({{ site.baseurl }}/_wiki/Structure-of-VBA's-generative-model.md), VBA handles a specific class of generative models, namely: non-linear state-space models. As we discuss below, this class of models is in fact very general, and can be used to emulate more sophisticated generative models that would apparently necessitate an extension of VBA's inference capabilities.


# Auto-regressive AR(1) models

First-order auto-regressive models or AR(1) capture the simplest form of stochastic processes. They are defined recusrvely as follows:

$$x_t= c+ \phi x_{t-1} +\epsilon_t$$

where $$c$$ and $$\phi$$ are constant parameters and $$\epsilon_t$$ is typically assumed to be i.i.d. Gaussian random noise.

> This process is [wide-sense stationary](https://en.wikipedia.org/wiki/Stationary_process#Weak_or_wide-sense_stationarity) if $$\mid \phi \mid <1$$.





# Auto-regressive AR(p) models

bla.

# Ornstein-Uhlenbeck processes

bla.

# GARCH models

bla.



