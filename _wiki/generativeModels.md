---
title: "How general is the class of VBA's generative models?"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

As we have described [here]({{ site.baseurl }}/_wiki/Structure-of-VBA's-generative-model.md), VBA handles a specific class of generative models, namely: non-linear state-space models. As we discuss below, this class of models is in fact very general, and can be used to emulate more sophisticated generative models that would apparently necessitate an extension of VBA's inference capabilities.


# Auto-regressive AR(1) models

First-order auto-regressive models or AR(1) capture the simplest form of stochastic processes. They are defined recusrvely as follows:

$$x_t= c+ \phi x_{t-1} +\eta_t$$

where $$c$$ and $$\phi$$ are constant parameters and $$\eta_t$$ is typically assumed to be i.i.d. Gaussian random noise.

> Some parameter constraints are necessary for the model to remain [wide-sense stationary](https://en.wikipedia.org/wiki/Stationary_process#Weak_or_wide-sense_stationarity). In particular, AR(1) processes with $$\mid \phi \mid >1$$ are not stationary.

Note that, in the seminal definition of AR(1) processes, the variable $$x_t$$ is directly observable. This corresponds to a specific case of nonlinear state-space models, where the evolution function of hidden states is given by the above equation, and the observation function is trivial and given by $$y_t = x_t + \epsilon_t$$ with a measure measurement noise precision $$\sigma \rightarrow \infty$$. Practically speaking, one can use high values for measure measurement noise precision (e.g. $$10^4$$).


# Ornstein-Uhlenbeck processes


bla.


# Auto-regressive AR(p) models

bla.



# GARCH models

bla.



