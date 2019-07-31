---
title: "Controlling the inversion options"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Strictly speaking, the main VBA inversion routine (`VBA_NLStateSpaceModel.m`) only requires the specification of evolution/observation functions and the model dimension (see [this page]({{ site.baseurl }}/wiki/VBA-model-inversion-in-4-steps) for a quick introduction).
However, optional arguments can be passed to this function through VBA's `options` input, which enables users to control the inversion. These are reviewed here.

# Dealing with categorical data

By default, VBA's generative model assumes observed data are [continuous](https://en.wikipedia.org/wiki/Continuous_function), for which there is a natural [distance metric](https://en.wikipedia.org/wiki/Metric_(mathematics)). Now if the data is categorical, there is no such natural metric, and one has to resort to probability distributions dealing with discrete events. For example, [binary](https://en.wikipedia.org/wiki/Binary_number) data can be treated as [binomial](https://en.wikipedia.org/wiki/Binomial_distribution) (Bernouilli) samples, whose sufficient statistic (first-order moment) is given by the observation function. This can be done by setting:

```matlab
options.sources.type = 1 ;
```

Note that this renders the measurement noise precision (and associated covariance components, see below) irrelevant. It turns out that this does not induce any major change in the VB inversion scheme under the Laplace approximation. In fact, when the dimension of the data is high enough, the empirical distribution of the "residuals" $$\epsilon = y - g(\vartheta)$$ will tend to a Gaussian density (this is actually used during VBA's initialization).

Data can also be categorical without being binary. In particular, it can have more than two "levels" (as for, e.g., colours or emotions, etc...). From a statistical perspective, this type of data can be modeled in terms of [multinomial](https://en.wikipedia.org/wiki/Multinomial_distribution) variables, which VBA can handle using:

```matlab
options.sources.type = 2 ;
```

In fact, the toolbox can concurrently fit a combination of multiple observations "sources" with different distributions and noise precisions. See [this page]({{site.baseurl }}/wiki/Multisources) for details regarding how to set up such a composite likelihood.

# Dealing with missing data

One may want to simply exclude data samples (either because one knows that they are unreliable, or because one has missing data). This can be done using the variable `options.isYout`. By default, `options.isYout` is a zero-valued matrix, whose size is identical to the data matrix `y`. Setting any of its value to one effectively asks VBA to disregard the corresponding data sample. For example:

```matlab
options.isYout(1:2,3) = 1 ;
```
means that the third time sample of the two first dimensions of `y` will not be considered in the inversion.

> **TIP:** when dealing with continuous data, the same can be obtained by setting the corresponding measurement precision matrix elements to zero: `options.priors.iQy{3}(1:2,1:2) = 0` (recall that by default, VBA assumes identity measurement precision matrices). Note that, rather than excluding a data sample, one may prefer to tell VBA about data samples' poor reliability, which can be done by resetting the corresponding measurement precision matrix elements to small but nonzero values...


# "Micro-time" resolution

The evolution equation is discrete in time. If one directly uses this as an approximation to a continuous dynamical system, then the approximation’s accuracy strongly depends on the data [sampling rate](https://en.wikipedia.org/wiki/Sampling_(signal_processing)#Sampling_rate). However, the toolbox allows one to specify a "microtime resolution", which is used to [recursively](https://en.wikipedia.org/wiki/Recursion_(computer_science)) apply the evolution function between two time samples. This can be useful to control the time [discretization errors](https://en.wikipedia.org/wiki/Discretization_error) introduced when approximating the original continuous dynamical system.

For example: how can one use a discrete time step of 100ms for the evolution function, given that the sampling rate is 1Hz?
Setting:

```matlab
options.decim = 10 ;
```
effectively asks VBA to recursively evaluate the evolution function ten times between two time samples. This resolves the problem if one passes a discrete time step of 100msec to the evolution function (see below).

One may also want to define the controlled input `u` time series at the micro-time resolution. This is done by setting:

```matlab
options.microU = 1 ;
```
If this is not the case, then VBA replicates the sampled input for each micro-time call of the evolution function.

# Optional arguments to evolution/observation functions

One may have to pass some arbitrary information to the evolution/observation functions. For example, when dealing with continuous dynamical systems, one may want to control the time discretization step (see above). This can be done by setting:

```matlab
options.inF.dt = 0.1 ;
```
where `inF.dt` is read out in the evolution function (and used accordingly; here: units would be seconds).

> **TIP:** Optional arguments to the observation function are passed through `options.inG`.

Note that the contents of the fields`options.inF` and `options.inG` are entirely arbitrary, since they are interpreted by the evolution/observation functions...

# Controlling the lagged Kalman forward pass

The VBA update of the hidden states is very similar in form to a [Kalman filter](https://en.wikipedia.org/wiki/Kalman_filter). More precisely, the scheme derives an approximation to the lagged posterior density $$p\left(x_t \mid y_{1:t+k}\right)$$, which contains the information about states at time $$t$$ given all observed data up to time $$t+k$$, where $$k$$ is the so-called **backward lag**. By default, this lag is set to 0. But it can be chosen arbitrarily, which allows one to infer on hidden states, whose changes impact observed data a few time samples later in time (e.g. due to some form of convolution operation). For example, the following sets the lag to 10:

```matlab
options.backwardLag = 10 ;
```

The main effect of increasing the lag is to average across more data points when deriving the hidden states, hence improving the precision (and the temporal smoothness) of the estimate.

> **TIP:** The ensuing [computational cost](https://en.wikipedia.org/wiki/Computational_complexity_of_mathematical_operations) scales with $$k^2 \times n^2$$, where $$n$$ is the states' dimension and $$k$$ is the Kalman backward lag. This means that the backward lag induces a trade-off between estimation efficiency and computational load. When setting the backward lag, we recommend to aim at finding the minimal value that still yields robust inference.


# Controlling VBA's algorithmic convergence

The free energy eventually serves as an analytical approximation to the log model evidence, but it is also used to monitor the algorithmic convergence of the algorithm. More precisely, VBA stops whenever the maximum number of iterations has been reached, or the increase in free energy has fallen below a predefined threshold. These two criteria can be controlled, e.g.:

```matlab
options.MaxIter = 10   ;
options.TolFun  = 1e-4 ;
```

Note that a minimum number of iterations can be fixed using `options.MinIter`.

In addition, VB updates of the posterior mean of hidden states and evolution/observation parameters rely upon an inner loop of iterative [Gauss-Newton](https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm) optimization schemes. Here again, the number of iterations as well as the relevant variational energy are monitored. They can be controlled by setting:

```matlab
options.GnMaxIter = 10   ;
options.GnTolFun  = 1e-4 ;
```

> **TIP:** instead of monitoring the variational energy through Gauss-Newton inner loops, one may want to force positive steps in the free energy itself. This can be done by setting:
>
>```matlab
options.gradF = 1 ;
```

# VBA initialization for stochastic dynamical systems

When dealing with [stochastic dynamical systems](http://www.scholarpedia.org/article/Stochastic_dynamical_systems), VBA initializes the inversion under a deterministic variant of the model. By default, this is done through an initial model inversion, where the prior mean of the states noise variance has been set to zero (with infinite prior precision, see [this page]({{ site.baseurl }}/wiki/Structure-of-VBA's-generative-model) for a generic description of VBA's priors). However, one may want to control the convergence of the VBA initialization. For example, setting:

```matlab
options.MaxIterInit = 0 ;
```
effectively bypasses the default VBA initialization.

In addition to initialized posterior densities on evolution/observation parameters and initial conditions, VBA initialization also derives the posterior $$p(\sigma\mid y,m)$$ on measurement noise precision. However, the inversion of the deterministic variant of the system might result in such strong model residuals, that one may want to ignore this. This can be done by setting:

```matlab
options.initHP = 0 ;
```
Here, the prior $$p(\sigma\mid m)$$ on measurement noise precision will be used instead of its posterior (under the deterministic variant model), when initializing the posterior for the main inversion of the stochastic system.

# Switching on/off graphical display

By default, VBA outputs a graphical summary of inversion results as the algorithm proceeds. However, this slows down the inversion. Setting:

```matlab
options.DisplayWin = 0 ;
```
effectively switches off the graphical display. NB: one can reproduce the graphical output post-hoc, by calling:

```matlab
VBA_ReDisplay(posterior,out) ;
```
where `posterior` and `out` are the two output arguments of the VBA inversion.

Having said this, one may want to eyeball inner Gauss-Newton iterations, for diagnosis purposes. This can be done by setting:

```matlab
options.GnFigs = 1 ;
```

> **TIP:** additional information is give in the main matlab window (e.g., default priors fill-in, free energy increases, etc...). This can also be switched off by setting:
>
> ```matlab
options.verbose = 0 ;
```


# Miscellaneous

## Bypassing the evolution function at times

One may want to skip the evolution function at times, i.e. replace it with the [identity mapping](https://en.wikipedia.org/wiki/Identity_function). For example, this is useful when dealing with learning models, where the transition from the initial conditions to the first states is ill-defined (because no feedback has yet been received). This can be done by setting the variable `options.skipf`. By default, `options.skipf` is a zero-valued matrix, whose length is the number of time samples. Setting any of its value to one effectively asks VBA to replace the corresponding states transition by the identity mapping. For example:

```matlab
options.skipf(1) = 1 ;
```
means that the first hidden states are copies of the initial conditions.

## Switching off some VB updates

One may want to bypass the VB updates of states' initial conditions and/or precision hyperparameters. This can be done by setting:

```matlab
options.updateX0 = 0 ;
```
and:

```matlab
options.updateHP = 0 ;
```
respectively.

> **TIP:** this is not identical to setting the variance of their respective priors to 0, because the free energy derivation still depends upon the prior and posterior uncertainty affecting these variables.

## Checking analytical gradients against numerical gradients

One may want to augment the evolution/observation functions with analytical gradients (w.r.t. states and parameters). This is useful to accelerate the inversion, because it bypasses the default numerical derivations. Because this may be tedious (and thus prone to errors), VBA allows model developpers to eyeball a comparison of their analytical gradients with numerical gradients. This is done by setting:

```matlab
options.checkGrads = 1 ;
```

## Dealing with delayed dynamical systems

VBA can deal with certain forms of [delays](https://en.wikipedia.org/wiki/Delay_differential_equation) by embedding the state space into an augmented state-space, which contains copies of states and their 'previous' values. The evolution function is then automatically wrapped, so that it receives delayed states. Such delays can be controlled using the variable `options.delays`. We refer the interested reader to the demonstration script `demo_delays.m`.
