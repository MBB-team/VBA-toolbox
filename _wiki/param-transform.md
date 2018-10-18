---
title: "Setting hard constraints through parameter transformations"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


In VBA, priors are Gaussian, i.e. parameters are allowed to vary without bounds (i.e. on the $$\left]-\infty,\infty \right[$$ interval). Setting informative priors allows one to insert "soft constraints" on parameter likely ranges. In particular, decreasing the prior variance eventually insures that parameter posterior estimates will be close to their prior mean. This is called "[shrinkage to the mean](https://en.wikipedia.org/wiki/Shrinkage_estimator)". This sort of constraint is "soft" because there always exist a certain amount of data that eventually allows posterior inference to escape from the influence of the priors. In other terms, as long as the prior variance is not zero, the strength of the prior constraint depends upon the informativeness of the data.

However, one may want to insert "[hard constraints](https://en.wikipedia.org/wiki/Constrained_optimization)", i.e. some models may include parameters that ought to be bounded within certain ranges. For example, [reinforcement learning](https://en.wikipedia.org/wiki/Reinforcement_learning) (RL) models are parameterized in terms of learning rates, which lie on the $$\left[0,1\right]$$ interval. Thus, to keep the parameters within bounds, one has to [re-parameterize](https://en.wikipedia.org/wiki/Parametrization) the model. Below are two overloaded examples of typical re-parameterizations for setting such hard constraints.



# Examples of re-parameterizations for setting hard constraints

## Range constraint

One may use the (re-scaled) [sigmoidal](https://en.wikipedia.org/wiki/Sigmoid_function) transform. For example, within a RL evolution function, the learning rate $$z$$ may be defined as $$z=s(x)$$, where $$x$$ is VBA's evolution parameter (which has a Gaussian prior). Thus, irrespective of the actual value of $$x$$, the learning rate $$z$$ will be constrained to lie on the $$\left[0,1\right]$$ interval.

## Positivity constraint

One may use the [exponential](https://en.wikipedia.org/wiki/Exponential_function) transform, i.e. $$exp(x)$$ will be positive, irrespective of the actual value of $$x$$. Another, less "stiff", positivity-enforcing mapping is $$exp(1+log(x))$$, which behaves linearly away from the origin.


> Note that Bayesian inference is **not invariant through changes in model parameterization**. For example, let us consider the mapping $$x \rightarrow x^3$$. This mapping does not insert any hard constraint, in the sense that the mapped parameter $$x^3$$ is allowed to vary without finite bounds. Nevertheless, Bayesian inference on a model that uses parameter $$x$$ does not yield the same result as inference on the same model, but this time using the mapped parameter $$x^3$$. This is because passing through the $$x \rightarrow x^3$$ can be seen as setting a non-Gaussian prior on $$x^3$$...


# Setting non-Gaussian priors

One may also want to impose a specific non-Gaussian prior to some parameters. For example, one may want to set [Beta](https://en.wikipedia.org/wiki/Beta_distribution) or [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) priors on model parameters. Although, superficially, VBA only deals with Gaussian priors, one can use parameter transformations to emulate non-Gaussian priors. This is also exemplified below.

## Setting a pre-specified prior density

Let $$z$$ be some model parameter, which one may want to set some prior density $$p_z(z)$$ on (e.g. Beta, Gamma, etc...). This can be done by mapping VBA's native parameter $$x$$ through the following transformation $$h(x)$$:

$$ h(x)=P_z\left(\Phi^{-1}\left(x\right)\right) $$

where $$\Phi(x)$$ is the cumulative distribution function of VBA's Gaussian priors on $$x$$ and $$P_z\left(z\right)$$ is the target cumulative distribution function of $$z$$ (ie. the integral of $$p_z(z)$$).

Note that the transformation $$h(x)$$ is the composition of two mappings:
- and the [inverse transform sampling](https://en.wikipedia.org/wiki/Inverse_transform_sampling) which produces a uniform distribution over [0,1] from any (here: gaussian) density,
- the [probability integral transform](https://en.wikipedia.org/wiki/Probability_integral_transform), which produces any distribution (here: $$p_z(z)$$) from the uniform distribution over [0,1].


> It turns ou that VBa's posterior inference does not depend *at all* on the way the native Gayssian prior is specified, as long as the tansformation $$h$$ uses the same normal cumulative distribution function...


## Emulating sparsity priors

So-called [sparse estimators](https://en.wikipedia.org/wiki/Sparse_approximation) arise in the context of model fitting, when one a priori assumes that only a few (unknown) model parameters deviate from zero. Sparsity constraints can be useful when the estimation problem is under-determined, i.e. when number of model parameters is much higher than the number of data points. In the context of [regularization](https://en.wikipedia.org/wiki/Regularization_(mathematics)) approaches, such constraints are enforced by minimizing [L1](http://mathworld.wolfram.com/L1-Norm.html) or [L0](https://en.wikipedia.org/wiki/Lp_space#When_p_=_0) norms, which yield the so-called [LASSO estimator](https://en.wikipedia.org/wiki/Lasso_(statistics)). 

Note that L1-norm minimization aprpoaches can be seen as a subcase of bayesian parameter estimation under non-Gaussian priors, more precisely: [Laplacian](https://en.wikipedia.org/wiki/Laplace_distribution) priors. In fact, this intuition generalizes to most sparsity constraints (cf. $$Lp$$-norm with $$p<2$$), which have their "sparsity prior" equivalent. It turns out that one can use the following simple parameter transform $$g_s$$ to emulate such sparsity priors without sacrificing the simplicity and robustness of VBA under Gaussian priors:

$$g_s(x)= \left(2 s(x) -1\right)x^2$$

where $$s$$ is the standard [sigmoid function](https://en.wikipedia.org/wiki/Sigmoid_function).

This "sparsify mapping" trick is demonstrated in this [technical note](https://arxiv.org/abs/1703.07168), which discloses the properties of such sparse estimators.


# Recovering the re-mapped posterior density

Finally, having performed a VBA analysis with some specified parameter transform in place, one may need to derive summary statistics for the transformed parameters (as opposed to directly using VBA's posterior pdf on native Gaussian parameters). This is the topic of the last section of this page...

## Laplace's method

If one is interested in the posterior density over transformed parameters, one can simply use the [Laplace approximation](https://en.wikipedia.org/wiki/Laplace's_method), which we summarize below.

Let $$f(x)$$ be the mapping used for setting a given hard constraint on some model parameter $$z=f(x)$$. Then a first-order [Taylor expansion](https://en.wikipedia.org/wiki/Taylor_series) yields:

$$ f(x) = f(E[x]) + f'(E[x])\times\left(\vartheta-E[x]\right) + ... $$

This first-order Taylor expansion is useful, because it can be used to derive the first- and second-order moments of $$z$$, given first- and second-order moments of $$x$$:

$$ E[z] \approx f(E[x]) $$

and

$$ V[z] \approx V[x]\times f'(E[x])^2 $$
 
The function `VBA_getLaplace` can be used to derive the above Laplace approximation, as follows:

```matlab
g_map = @myMapping;
dim = struct('n',0,'n_theta',0),'n_phi',1);
opt.priors = posterior;
[muP,VP] = VBA_getLaplace([],[],g_map,dim,opt,0)
```

where `@myMapping` implements the parameter tansformation (but with the usual VBA i/o) and `posterior` has been obtained using VBA...


## Monte-Carlo method

 Alternatively, one can [sample](https://en.wikipedia.org/wiki/Monte_Carlo_method) (see `VBA_sample.m`) from the Gaussian posterior density over un-transformed parameters, pass the samples through the transform, and then report [summary statistics](https://en.wikipedia.org/wiki/Summary_statistics) over the set of transformed samples (such as mean and variance) and/or full sampling histograms. Such sampling approach can be used, for example, in the aim of recovering [credible intervals](https://en.wikipedia.org/wiki/Credible_interval) over constrained (mapped) parameters.




