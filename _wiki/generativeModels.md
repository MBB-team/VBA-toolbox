---
title: "How general is the class of VBA's generative models?"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

As we have described [here]({{ site.baseurl }}/_wiki/Structure-of-VBA's-generative-model.md), VBA handles a specific class of generative models, namely: non-linear state-space models. As we discuss below, this class of models is in fact very general, and can be used to emulate more sophisticated generative models that would apparently necessitate an extension of VBA's inference capabilities.


# Auto-regressive AR(1) models

First-order auto-regressive models or AR(1) capture the simplest form of [markovian stochastic processes](https://en.wikipedia.org/wiki/Markov_chain), namely: [random walks](https://en.wikipedia.org/wiki/Random_walk). They are defined recursively as follows:

$$x_t= c+ \theta x_{t-1} +\eta_t$$

where $$c$$ and $$\theta$$ are constant parameters and $$\eta_t$$ is typically assumed to be i.i.d. Gaussian random noise.

> Some parameter constraints are necessary for the model to remain [wide-sense stationary](https://en.wikipedia.org/wiki/Stationary_process#Weak_or_wide-sense_stationarity). In particular, AR(1) processes with $$\mid \phi \mid >1$$ are not stationary.

Note that, in the seminal definition of AR(1) processes, the variable $$x_t$$ is directly observable. This corresponds to a specific case of nonlinear state-space models, where the evolution function of hidden states is given by $$f(x_t)=c+\theta x_t$$, and the observation function is trivial and given by $$g(x_t) = x_t$$ with a measure measurement noise precision $$\sigma \rightarrow \infty$$. Practically speaking, one can use high values for measure measurement noise precision (e.g. $$10^4$$).


# Ornstein-Uhlenbeck processes

[Ornstein-Uhlenbeck processes](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process) are a continuous variant of AR(1) processes, which tend to drift towards its long-term mean (which is why it called [mean-reverting](https://en.wikipedia.org/wiki/Mean_reversion_(finance))). Formally speaking, an Ornstein–Uhlenbeck process $$x(t)$$ satisfies the following [stochastic differential equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation):

$$dx(t) = \theta\left(\mu-x(t)\right) + \beta dW(t)$$

where $$\theta$$ and $$\beta$$ are constant parameters and $$dW(t)$$ is a (Wiener process](https://en.wikipedia.org/wiki/Wiener_process).

VBA can approximate such continuous process using the following evolution function: $$f(x_t)=x_t+\delta t \left(\mu-x_t\right)$$, where $$\delta t$$ is the discretization step. 







# Auto-regressive AR(p) models

bla.



# GARCH models

bla.



