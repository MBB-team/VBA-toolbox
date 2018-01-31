---
title: "Extending VBA's generative model"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

As we have described [here]({{ site.baseurl }}/wiki/Structure-of-VBA's-generative-model), VBA deals with so-called *nonlinear state-space models*. As we discuss below, this class of generative models is in fact very general, and can be used to emulate more sophisticated generative models that would apparently necessitate an extension of VBA's inference capabilities. What follows is a mathematical note that exemplifies how to frame apparently VBA-incompatible statistical models as state-space models that VBA can handle.


# Auto-regressive AR(1) models

First-order auto-regressive models or AR(1) capture the simplest form of [markovian stochastic processes](https://en.wikipedia.org/wiki/Markov_chain), namely: [random walks](https://en.wikipedia.org/wiki/Random_walk). They are defined recursively as follows:

$$x_t= c+ \theta x_{t-1} +\eta_t$$

where $$c$$ and $$\theta$$ are constant parameters, and $$\eta_t$$ is typically assumed to be i.i.d. Gaussian random noise.

> Some parameter constraints are necessary for the model to remain [wide-sense stationary](https://en.wikipedia.org/wiki/Stationary_process#Weak_or_wide-sense_stationarity). In particular, AR(1) processes with $$\mid \theta \mid >1$$ are not stationary.

Note that, in the seminal definition of AR(1) processes, the variable $$x_t$$ is directly observable. This corresponds to a specific case of nonlinear state-space models, where the evolution function of hidden states is given by $$f(x_t)=c+\theta x_t$$, and the observation function is trivial and given by $$g(x_t) = x_t$$ with a measurement noise precision $$\sigma \rightarrow \infty$$. Practically speaking, one can use high values for measurement noise precision (e.g. $$10^4$$).


# Ornstein-Uhlenbeck processes

An [Ornstein-Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process) is a continuous variant of AR(1) processes that drifts back and forth towards its long-term mean, which is why it is called [mean-reverting](https://en.wikipedia.org/wiki/Mean_reversion_(finance)). Formally speaking, an Ornstein–Uhlenbeck process $$x(t)$$ satisfies the following [stochastic differential equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation):

$$dx(t) = \theta\left(\mu-x(t)\right) + \beta dW(t)$$

where $$\theta$$ and $$\beta$$ are constant parameters, $$\mu$$ is the long-term mean of the process, and $$dW(t)$$ is a [Wiener process](https://en.wikipedia.org/wiki/Wiener_process). As can be seen in the equation above, the process is expected to exhibit some form of [regression to the mean](https://en.wikipedia.org/wiki/Regression_toward_the_mean), because deviations from the mean $$\mu-x(t)$$ effectively induce restoring forces.

VBA can approximate such continuous process using the following evolution function: $$f(x_t)=x_t+\Delta t \theta\left(\mu-x_t\right)$$, where $$\Delta t$$ is the discretization step. Evolution parameters now include both $$\mu$$ and $$\theta$$. The continuous limit is obtained by increasing the number of recursive calls to the evolution function between two time samples (e.g., by setting VBA's micro-time resolution to `options.decim = 10` or more). Note that the variance $$\beta$$ of the approximated Wiener process can be recovered from VBA's estimate of the state noise precision.


# Auto-regressive AR(p) models

Higher-order auto-regressive models or [AR(p)](https://en.wikipedia.org/wiki/Autoregressive_model) are a non-Markovian generalization of AR(1) processes:

$$ x_t= c+ \sum_{i=1}^{p}\theta_i x_{t-i} +\eta_t $$

where $$c$$ and $$\theta_1,...,\theta_p$$ are constant parameters, and $$\eta_t$$ is assumed to be i.i.d. Gaussian random noise.

In principle, VBA only deals with [Markovian systems](https://en.wikipedia.org/wiki/Markov_chain), i.e. systems whose evolution depends solely upon their current state. Thus, the non-Markovian structure of AR(p) processes should eschew any VBA-based system identification. But in fact, a very simple solution to this apparent issue is to augment the native state-space of AR(p) processes with dummy states that are copies of past instances of native states. As we will see, this allows us to describe any non-Markovian system in terms of a (higher-dimensional) Markovian system. Let us define $$z_t$$ as the following $$p\times 1$$ augmented state variable:

$$ z_t = \left[\begin{array}{l} x_t \\ x_{t-1} \\ \vdots \\ x_{t-p+1} \end{array}\right] \implies z_{t+1} = \left[\begin{array}{l} x_{t+1} \\ x_{t} \\ \vdots \\ x_{t-p+2}\end{array}\right] = f(z_t) $$

where $$p$$ is the target order of the autoregressive process. This augmented state-space can be used to circumvent non-Markovian dynamics of AR(p) processes, the structure of which can be emulated using the following evolution function on $$z_t$$:

$$ f(z_t) = \left[\begin{array}{l} c+ \sum_{i=1}^{p}\theta_i x_{t-i+1} \\ x_{t} \\ \vdots \\ x_{t-p+2} \end{array}\right] = \left[\begin{array}{l} c \\ 0 \\ \vdots \\ 0 \end{array}\right] + \left[\begin{array}{l} A^T \\ {L_1}^T \\ \vdots \\ {L_{p-1}}^T \end{array}\right] z_t  $$

where $$A$$ and $$L_1,...,L_{p-1}$$ are $$p\times 1$$ vectors given by:

$$ A = \left[\begin{array}{l} \theta_1 \\ \theta_2 \\ \vdots \\ \theta_p \end{array}\right], L_1 = \left[\begin{array}{l} 1 \\ 0 \\ \vdots \\ 0 \end{array}\right], ..., L_{p-1} = \left[\begin{array}{l} 0 \\ \vdots \\ 1 \\ 0 \end{array}\right] $$

One can see that the native unidimensional AR(p) process is now described in terms of Markovian dynamics on a p-dimensional state-space $$z$$, whose evolution function $$f(z)$$ is compatible with VBA analyses.

The corresponding observation function would then simply be given by:

$$ g(z_t) = {L_1}^T z_t = x_t$$

with a measurement noise precision $$\sigma \rightarrow \infty$$  (in practice, $$10^4$$ or so).

> Evolution parameters include $$c$$ and $$\theta_1,...,\theta_p$$, which can be estimated using VBA in the usual way. Note that in this example, there is no observation parameter.

Now, there is a last problem to fix. Recall that, when inverting stochastic models, VBA assumes that hidden states's dynamics are driven by a mixture of deterministic (the evolution function) and ramdom forces (state noise), i.e.: $$z_{t+1} = f(z_t) + \eta_t$$. State noise is required for AR(p) processes because it eventually triggers observed variations in $$x_t$$. However, we do not want state noise to perturb the "copy-paste" operation performed on the $$p-1$$ last entries of $$z_t$$. In principle, this can be achieved by resetting the state noise precision matrix $${Q_x}^{-1}$$ as the following diagonal matrix:

$$ {Q_x}^{-1} = \left[\begin{array}{cccc} 1 & 0 & \cdots & 0 \\ 0 & r & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & r  \end{array}\right] $$

where $$r \rightarrow \infty $$ to ensure that the past history of hidden states is transfered without distortion. Practically speaking, this can be approximated by changing the matrix `options.priors.iQx{t}` as above, with $$r$$ set to an appriately high value (e.g., $$10^4$$).

Finally, we recommend that the *[VBA's Kalman backward lag]({{ site.baseurl }}/wiki/Controlling-the-inversion-using-VBA-options/#controlling-the-lagged-kalman-forward-pass)* be adjusted to the order of the AR(p) process, as follows: `options.backwardLag=p`. This is necessary to properly account for the delayed influence of the states' history.

> It is straightforward to generalize the above "augmented state-space" trick to multivariate AR(p) proceses, or, in fact, to any form of non-Markovian dynamics (including, e.g., delayed dynamical systems...)!


# ARCH models

In brief, [autoregressive conditional heteroskedastic (ARCH)](https://en.wikipedia.org/wiki/Autoregressive_conditional_heteroskedasticity#GARCH) processes are such that the variance of the current error term or innovation depends upon the system's state. In discrete time, an example of a generalized ARCH process would be given by:

$$x_t= a\left(x_{t-1}\right) + \beta\left(x_{t-1}\right)\eta_t$$

where $$a\left(x\right)$$ is an arbitrary mapping of system's states and $$\beta\left(x\right)$$ acts as the state-dependent standard deviation of state noise $$\eta$$. Note that this dependency can be arbitrary, which endows ARCH models with a great deal of flexibility. Typically, ARCH models are employed in modeling financial time series that exhibit time-varying [volatility](https://en.wikipedia.org/wiki/Stochastic_volatility), i.e. periods of swings (when $$\beta$$ is high) interspersed with periods of relative calm (when $$\beta$$ is low).

In its native form, VBA's generative model does not apprently handle state-dependent noise. But in fact, one can use the same trick as for AR(p) models. First, one augments the native state-space with dummy states $$w_t$$ as follows:

$$ z_t = \left[\begin{array}{l} x_t \\ w_t \end{array}\right] $$

As we will see, the dummy state $$w_t$$ will be composed of "pure" noise.

Second, one defines their evolution function as follows:

$$ f(z_t)  = \left[\begin{array}{l} a\left(x_t\right) + \beta\left(x_t\right)w_t \\ 0 \end{array}\right] = \left[\begin{array}{l} a\left({L_1}^T z_t\right) + \beta\left({L_1}^T z_t\right) {L_2}^T z_t \\ 0 \end{array}\right] $$

where $$L_1$$ and $$L_2$$ are are $$2\times 1$$ vectors given by:

$$ L_1 = \left[\begin{array}{l} 1 \\ 0 \end{array}\right],L_2 = \left[\begin{array}{l} 0 \\ 1 \end{array}\right] $$

Since the deterministic flow of dummy states is null, their stochastic dynamics are solely driven by state noise, i.e. $$w_t = \eta_t$$. 

> Both the native deterministic flow $$a(x)$$ and the state-dependent standard-deviation $$\beta(x)$$ can be parameterized through evolution parameters, i.e.: $$a(x)=a(x,\theta)$$ and $$\beta(x)=\beta(x,\theta)$$. These can then be estimated using VBA.

Third, one resets the state noise precision matrix $${Q_x}^{-1}$$ as follows:

$$ {Q_x}^{-1} = \left[\begin{array}{cc} r & 0 \\ 0 & 1  \end{array}\right] $$

with $$r \rightarrow \infty $$. This ensures that state noise $$\eta$$ only enters at the level of dummy states $$w$$, which is then rescaled by $$\beta\left(x\right)$$ before perturbing hidden states. Practically speaking, this can be approximated by changing the matrix `options.priors.iQx{t}` as above, with $$r$$ set to an appriately high value (e.g., $$10^4$$).

Finally, one defines the observation function on the augmented state-space as follows:

$$ g(z_t) = {L_1}^T z_t = x_t$$

with a measure measurement noise precision $$\sigma \rightarrow \infty$$ (in practice, $$10^4$$ or so).

Note that, as for AR(p) processes, we recommend that the *[VBA's Kalman backward lag]({{ site.baseurl }}/wiki/Controlling-the-inversion-using-VBA-options/#controlling-the-lagged-kalman-forward-pass)* be increased as much as possible (at least `options.backwarLag=2`, but preferably more...).

> The computational cost incurred when emulating ARCH models using homoscedastic state-space models is twofold: (i) some form of nonlinearity in the evolution function, and (ii) an increase in the dimensionality of the system. Note that the latter cost is somehow paid *twice* here, since the computational load of VBA's inversion scales with the Kalman backward lag.


# Covariance component models

For the sake of simplicity, we will only consider below static generative models of the form $$y=g(\phi)+\eta$$, where $$y$$ is the data, $$ \phi$$ are uknown observation parameters, $$g$$ is an arbitrary observation function and $$\eta$$ are model residuals.

Recall that, by default, VBA can only estimate one variance hyperparameter (namely: $$\sigma$$) per data source. In other terms, the covariance of model residuals $$\eta$$ is constrained to be a rescaling of a fixed covariance matrix $$Q_y$$, i.e.: $$E[\eta\eta^T]= \sigma^{-1}Q_y$$. However, one may have prior information regarding the statistical structure of model residuals, in the form of a mixture of covariance components:

$$E[\eta\eta^T]= \sum_i \lambda_i Q_i$$

where $$Q_i$$ are fixed covariance components and $$\lambda_i$$ are unknown mixture coefficients. Covariance component models simply provide estimates of both model parameters $$\phi$$ and hyperparameters $$\lambda$$. Critical here is the fact that structured residuals can (and should) change parameters estimates...

Here as well, one can use a similar trick as above. To begin with, note that the above covariance component mixture would follow from defining model residuals $$\eta$$ as the following weighted sum of dummy random variables $$w_i$$:

$$ \Bigg\{ \begin{array}{l} \eta = \sum_i \sqrt{\lambda_i} w_i \\ E[w_i w_i^T]= Q_i \end{array} $$

Thus, it suffices to augment the native parameter space with two sets of dummy parameters $$z$$ and $$\lambda$$, and replace the native observation function $$g$$ by the following augmented observation function $$h$$:

$$ h(\phi,z,\lambda) = g(\phi) + \sum_i \sqrt{\lambda_i} U_i z_i$$

where $$U_i$$ are the known [matricial square root](https://en.wikipedia.org/wiki/Square_root_of_a_matrix) of covariance components $$Q_i$$, i.e.: $$Q_i=U_i U_i^T$$ (these can be obtained from numerical [SVD decompositions](https://en.wikipedia.org/wiki/Singular_value_decomposition)).

Setting i.i.d. Gaussian priors on dummy variables $$z$$ would then emulate covariance component models. On a practical note, parameter estimation would be better behaved if one used some form of [hard positivity constraint]({{ site.baseurl }}/wiki/param-transform) on the $$\lambda$$s.

> If the covariance components reduce to channel-specific variances, then one can use VBA's ["multi-source" inversion]({{ site.baseurl }}/wiki/Multisources) as a simple and elegant shortcut!

Note that, by construction, native model parameters $$\phi$$ and dummy noise variables $$z$$ compete for explaining variability in observed data $$y$$. This is not an artefactucal consequence of our way of treating covariance component models. This competition is simply more implicit in the usual model inversion framework, whereby one does not directly derive posterior densities over noise variables, but rather provide estimates from the prediction error $$y-g(\phi)$$...



# Switching dynamical systems

Many natural systems, ranging from neurons firing patterns to collective motion of animal crowds, give rise to time series data with complex dynamics that may not be amenable to explicit modelling. Neverthless, one can gain insight into these systems by decomposing the data into segments that are each explained by simpler dynamic units. This induces so-called switching dynamical systems. The state-space of such system is composed of both continuous ($$x$$) and discrete ($$z$$) states. The latter control the form of the evolution function of the former, eventually capturing potentially abrupt changes in the system's deterministic flow.

Strictly speaking, a switching dynamical system cannot be derived as a limiting subcase of nonlinear state-space models. The issue here, derives from (mutually exclusive) discrete states, that evolve according to transition probabilities that have no formal equivalent in continuous state-spaces...

Having said this, one can approximate discrete states in terms of continuous states passed through steep sigmoidal mappings. For example, let $$z$$ be the $$n\times 1$$ [indicator vector](https://en.wikipedia.org/wiki/Indicator_vector) of a discrete state, i.e. the only non-zero entry of $$z$$ specifies which element of the discrete space is currently active. Then, there always exists some continuous $$n\times 1$$ vector $$x$$ that verifies:

$$z = s_{\beta}\left(x\right) = \frac{\exp \beta x}{\sum_{i=1}^n \exp \beta x_i}$$

where $$s_{\beta}()$$ is a softmax mapping with inverse-temperature $$\beta \rightarrow \infty $$.

In switching dynamical systems, such discrete states evolve according to the following transition probability matrix:

$$ P\left(z_t\mid z_{t-1}\right) = \Pi = \left[\begin{array}{cccc} \pi_{11} & \pi_{12} & \cdots & \pi_{1n} \\ \pi_{21} & \pi_{22} & \cdots & \pi_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ \pi_{n1} & \pi_{n2} & \cdots & \pi_{nn} \end{array}\right] $$

where $$1 = \sum_{i=1}^{n} \pi_{ij}$$ for all $$j$$.

Conditional on $$z_{t-1}$$, the first- and second-order moments of $$z_t$$ are thus given by:

$$ E\left[z_t\mid z_{t-1}\right] = \Pi z_{t-1} $$

and:

$$ V\left[z_t\mid z_{t-1}\right] = \textsf{diag}(\Pi z_{t-1}) - \Pi z_{t-1} \Pi z_{t-1}^T $$

where $$\textsf{diag}(x)$$ is, by abuse of notation, the $$n\times n$$ matrix whose diagonal entries is composed of the vector $$x$$.

Thus, switching dynamical systems can be seen as some form of multivariate ARCH model, whereby the noise variance-covariance matrix is dependent upon the previous state of the system (and given by the above equation). One would then define an augmented state-space using $$n$$ continuous states mapped through a sigmoid mapping and $$n$$ dummy states that are driven by noise but are rescaled approrpiately prior to perturbing the first half of the state-space...





