---
title: "Extending VBA's generative models"
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

VBA can approximate such continuous process using the following evolution function: $$f(x_t)=x_t+\Delta t \left(\mu-x_t\right)$$, where $$\Delta t$$ is the discretization step. The continuous limit is obtained by increasing the number of recursive calls to the evolution function between two time samples (e.g., by setting VBA's micro-time resolution to `options.decim = 10`). 


# Auto-regressive AR(p) models

Higher-order auto-regressive models or [AR(p)](https://en.wikipedia.org/wiki/Autoregressive_model) are a non-Markovian generalization of AR(1) processes:

$$ x_t= c+ \sum_{i=1}^{p}\theta_i x_{t-i} +\eta_t $$

where $$c$$ and $$\theta_1,...,\theta_p$$ are constant parameters and $$\eta_t$$ is assumed to be i.i.d. Gaussian random noise.

The only issue here stems from the non-Markovian structure of AR(p) processes. But in fact, a very simple solution to this is to augment the state-space with dummy states that are copies of past instances of native states. Let us define $$z_t$$ as the following dummay state variable:

$$ z_t = \left[\begin{array}{l} x_t \\ x_{t-1} \\ \vdots \\ x_{t-p+1} \end{array}\right] \implies z_{t+1} = \left[\begin{array}{l} x_{t+1} \\ x_{t} \\ \vdots \\ x_{t-p+2}\end{array}\right] $$

where $$p$$ is the target order of the autoregressive process. Then the structure of AR(p) processes can be emulated using the following evolution function on $$z_t$$:

$$ f(z_t) = \left[\begin{array}{l} c \\ 0 \\ \vdots \\ 0 \end{array}\right] + \left[\begin{array}{l} A^T \\ {L_1}^T \\ \vdots \\ {L_{p-1}}^T \end{array}\right] z_t = \left[\begin{array}{l} c+ \sum_{i=1}^{p}\theta_i x_{t-i+1} \\ x_{t} \\ \vdots \\ x_{t-p+2} \end{array}\right] $$

where $$A$$ and $$L_1,...,L_{p-1}$$ are $$p\times 1$$ vectors given by:


$$ A = \left[\begin{array}{l} \theta_1 \\ \theta_2 \\ \vdots \\ \theta_p \end{array}\right], L_1 = \left[\begin{array}{l} 1 \\ 0 \\ \vdots \\ 0 \end{array}\right], ..., L_{p-1} = \left[\begin{array}{l} 0 \\ \vdots \\ 0 \\ 1 \end{array}\right] $$

The corresponding observation function would then be simply given by:

$$ g(z_t) = {L_1}^T z_t = x_t$$

with a measure measurement noise precision $$\sigma \rightarrow \infty$$.

Now, there is a last problem to fix. Recall that, when inverting stochastic models, VBA assumes that hidden states's evolution is driven by a mixture of deterministic (the evolution function) and ramdom forces (state noise), i.e.: $$z_{t+1} = f(z_t) + \eta_t$$. State noise is required for AR(p) processes because it eventually triggers observed variations in $x_t$. However, we do not want state noise to perturb the "copy-paste" operation performed on the $p-1$ last entries of $z_t$. In principle, this can be achieved by resetting the state noise precision matrix $${Q_x}^{-1}$$ as the following diagonal matrix:

$$ {Q_x}^{-1} = \left[\begin{array}{cccc} 1 & 0 & \cdots & 0 \\ 0 & r & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & r  \end{array}\right] $$

where $$r \rightarrow \infty $$. Practically speaking, this can be approximated by changing the matrix `options.priors.iQx{t}` as above, with $$r$$ set to an appriately high value (e.g., $$10^4$$).


# ARCH models

In brief, [autoregressive conditional heteroskedastic (ARCH)](https://en.wikipedia.org/wiki/Autoregressive_conditional_heteroskedasticity#GARCH) processes are such that the variance of the current error term or innovation depends upon the system's state. In discrete time, an example of a generalized ARCH process would be given by:

$$x_t= a\left(x_{t-1}\right) + \beta\left(x_{t-1}\right)\eta_t$$

where $$a\left(x\right)$$ is an arbitrary mapping of system's states and $$\beta\left(x\right)$$ acts as the state-dependent standard deviation of state noise. Note that this dependency can be arbitrary, which endows ARCH models with a great deal of flexibility. Typically, ARCH models are employed in modeling financial time series that exhibit time-varying [volatility](https://en.wikipedia.org/wiki/Stochastic_volatility), i.e. periods of swings (when $$\beta$$ is high) interspersed with periods of relative calm (when $$\beta$$ is low).

In its native form, VBA's generative model does not apprently handle state-dependent noise. But in fact, one can use the same trick as for AR(p) models. In brief:


First, one augments hidden states with dummy states $$w_t$$ as follows:

$$ z_t = \left[\begin{array}{l} x_t \\ w_t \end{array}\right] $$

As we will see, the dummy state $$w_t$$ will be composed of "pure" noise.

Second, one defines their evolution function as follows:

$$ f(z_t) = \left[\begin{array}{l} a\left({L_1}^T z_t\right) + \beta\left({L_1}^T z_t\right) {L_2}^T z_t \\ 0 \end{array}\right] = \left[\begin{array}{l} a\left(x_{t-1}\right) + \beta\left(x_{t-1}\right)w_t \\ 0 \end{array}\right] $$

where $$L_1$$ and $$L_2$$ are are $$2\times 1$$ vectors given by:

$$ L_1 = \left[\begin{array}{l} 1 \\ 0 \end{array}\right],L_2 = \left[\begin{array}{l} 0 \\ 1 \end{array}\right] $$

In turn, VBA assumes that dummy states are entirely driven by state noise, i.e. $$w_t = \eta_t$$.

Third, one resets the state noise precision matrix $${Q_x}^{-1}$$ as follows:

$$ {Q_x}^{-1} = \left[\begin{array}{cc} r & 0 \\ 0 & 1  \end{array}\right] $$

with $$r \rightarrow \infty $$. This ensures that state noise $$\eta$$ only enters at the level of dummy states $$w$$, which is then rescaled by $$\sigma\left(x\right)$$ before perturbing hidden states $$x$$. Practically speaking, this can be approximated by changing the matrix `options.priors.iQx{t}` as above, with $$r$$ set to an appriately high value (e.g., $$10^4$$).

Finally, one defines the observation function on the augmented state-space as follows:

$$ g(z_t) = {L_1}^T z_t = x_t$$

with a measure measurement noise precision $$\sigma \rightarrow \infty$$.


# Switching dynamical systems

Many natural systems, ranging from neurons firing patterns to collective motion of animal crowds, give rise to time series data with complex, nonlinear dynamics. One can gain insight into these systems by decomposing the data into segments that are each explained by simpler dynamic units. This induces so-called switching dynamical systems. The state-space of such system is composed of both continuous ($$x$$) and discrete ($$z$$) states. In particular, the discrete states control the form of the evolution function of the continuous states, eventually inducing abrupt changes in the system's deterministic flow.

Strictly speaking, a switching dynamical system cannot be derived as a limiting subcase of nonlinear state-space models. The issue here, derives from (mutually exclusive) discrete states, that evolve according to transition probabilities that have no formal equivalent in continuous state-spaces...

Having said this, one can approximate discrete states in terms of continuous states passed through steep sigmoidal mappings. For example, let $$z$$ be the $$n\times 1$$ [indicator vector](https://en.wikipedia.org/wiki/Indicator_vector) of a discrete state, i.e. the only non-zero entry of $$z$$ specifies which element of the discrete space is currently active. Then, there always exists some continuous $$n\times 1$$ vector $$x$$ that verifies:

$$z = s_{\beta}\left(x\right) = \frac{\exp \beta x}{\sum_{i=1}^n \exp \beta x_i}$$

where $$s_{\beta}()$$ is a softmax mapping with inverse-temperature $$\beta \rightarrow \infty $$.

In switching dynamical systems, such discrete states evolve according to the following transition probability matrix:

$$ P\left(z_t\mid z_{t-1}\right) = \Pi = \left[\begin{array}{cccc} \pi_{11} & \pi_{12} & \cdots & \pi_{1n} \\ \pi_{21} & \pi_{22} & \cdots & \pi_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ \pi_{n1} & \pi_{n2} & \cdots & \pi_{nn} \end{array}\right] $$

where $$1 = \sum_{i=1}^{n} \pi_{ij}$$ for all $$j$$.

Conditional on $$z_{t-1}$$, the first-order and second-order moment of $$z_t$$ are thus given by:

$$ E\left[z_t\mid z_{t-1}\right] = \Pi z_{t-1} $$

and:

$$ V\left[z_t\mid z_{t-1}^{(j)}=1\right] = diag(\Pi z_{t-1}) - \left(\Pi z_{t-1}\right) \left(\Pi z_{t-1}\right)^T $$

where $$diag(x)$$ is, by abuse of notation, the $$n\times n$$ matrix whose diagonal entries is composed of the vector $$x$$.

Thus, switching dynamical systems can be seen as some form of multivariate ARCH model, whereby the noise variance-covariance matrix is dependent upon the previous state of the system (and given by the above equation). One would then define a state-space using $$n$$ continuous states mapped through a sigmoid mapping and $$n$$ dummy states that are driven by noise but are rescaled prior to perturbing the first half of the state-space...





