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

An [Ornstein-Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process) is a continuous variant of AR(1) processes that tends to drift towards its long-term mean, which is why it is called [mean-reverting](https://en.wikipedia.org/wiki/Mean_reversion_(finance)). Formally speaking, an Ornstein–Uhlenbeck process $$x(t)$$ satisfies the following [stochastic differential equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation):

$$dx(t) = \theta\left(\mu-x(t)\right) + \beta dW(t)$$

where $$\theta$$ and $$\beta$$ are constant parameters, $$\mu$$ is the long-term mean of the process, and $$dW(t)$$ is a [Wiener process](https://en.wikipedia.org/wiki/Wiener_process). As can be seen in the equation above, the process is expected to exhibit some form of [regression to the mean](https://en.wikipedia.org/wiki/Regression_toward_the_mean), because deviations from the mean $$\mu-x(t)$$ effectively induce restoring forces.

VBA can approximate such continuous process using the following evolution function: $$f(x_t)=x_t+\delta t \left(\mu-x_t\right)$$, where $$\delta t$$ is the discretization step. The continuous limit is obtained by increasing the number of recursive calls to the evolution function between two time samples (e.g., by setting VBA's micro-time resolution to `options.decim = 10`). 


# Auto-regressive AR(p) models

Higher-order auto-regressive models or [AR(p)](https://en.wikipedia.org/wiki/Autoregressive_model) are a non-Markovian generalization of AR(1) processes:

$$x_t= c+ \sum_{i=1}^{p}\theta_i x_{t-i} +\eta_t$$

where $$c$$ and $$\theta_1,...,\theta_p$$ are constant parameters and $$\eta_t$$ is assumed to be i.i.d. Gaussian random noise.

The only issue here stems from the non-Markovian structure of AR(p) processes. But in fact, a very simple solution to this is to augment the state-space with dummy states that are copies of past instances of native states. Let us define $$z_t$$ as the following dummay state variable:

$$ z_t = \left[\begin{array}{l} x_t \\ x_{t-1} \\ ...\\ x_{t-p+1} \end{array}\right]$$

where $$p$$ is the target order of the autoregressive process. Then the structure of AR(p) processes can be emulated using the following evolution function on $$z_t$$:

$$ f(z_t) = \left[\begin{array}{l} c \\ 0 \\ ...\\ 0 \end{array}\right] + \left[\begin{array}{l} A^T \\ {L_1}^T \\ ...\\ {L_{p-1}}^T \end{array}\right] z_t $$

where $$A$$ and $$L_1,...,L_{p-1}$$ are vectors given by:

$$ A = \left[\begin{array}{l} \theta_1 \\ \theta_2 \\ ...\\ \theta_p \end{array}\right], L_1 = \left[\begin{array}{l} 0 \\ 1 \\ ...\\ 0 \end{array}\right], ..., L_{p-1} = \left[\begin{array}{l} 0 \\ 0 \\ ...\\ 1 \end{array}\right]$$.




# GARCH models

bla.



