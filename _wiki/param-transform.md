---
title: "Setting hard constraints through parameter transformations"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


In VBA, priors are Gaussian, i.e. parameters are allowed to vary without bounds (i.e. on the $$\left]-\infty,\infty \right[$$ interval). Setting informative priors allows one to insert "soft constraints" on parameter likely ranges. In particular, decreasing the prior variance eventually insures that parameter posterior estimates will be close to their prior mean. This is called "[shrinkage to the mean](https://en.wikipedia.org/wiki/Shrinkage_estimator)". This sort of constraint is "soft" because there always exist a certain amount of data that eventually allows posterior inference to escape from the influence of the priors. In other terms, as long as the prior variance is not zero, the strength of the prior constraint depends upon the informativeness of the data.

However, one may want to insert "[hard constraints](https://en.wikipedia.org/wiki/Constrained_optimization)", i.e. some models may include parameters that ought to be bounded within certain ranges. For example, [reinforcement learning](https://en.wikipedia.org/wiki/Reinforcement_learning) (RL) models are parameterized in terms of learning rates, which lie on the $$\left[0,1\right]$$ interval. Thus, to keep the parameters within bounds, one has to [re-parameterize](https://en.wikipedia.org/wiki/Parametrization) the model. Below are two overloaded examples of typical re-parameterizations for setting such hard constraints.



# Examples of re-parameterizations for setting hard constraints

## range constraint

One may use the (re-scaled) [sigmoidal](https://en.wikipedia.org/wiki/Sigmoid_function) transform. For example, within a RL evolution function, the learning rate may be defined as `sigm(P)`, where `P` is VBA's evolution parameter (which has a Gaussian prior). Thus, irrespective of the actual value of `P`, the learning rate `sigm(P)` will be constrained to lie on the $$\left[0,1\right]$$ interval.

## Positivity constraint

One may use the [exponential](https://en.wikipedia.org/wiki/Exponential_function) transform, i.e. `exp(P)` will be positive, irrespective of the actual value of `P`. Another, less "stiff", positivity-enforcing mapping is `exp(1+log(x))`, which behaves linearly away from the origin.


> Note that Bayesian inference is **not invariant through changes in model parameterization**. For example, let us consider the mapping $$P \rightarrow P^3$$. This mapping does not insert any hard constraint, in the sense that the mapped parameter $$P^3$$ is allowed to vary without finite bounds. Nevertheless, Bayesian inference on a model that uses parameter `P` does not yield the same result as inference on the same model, but this time using the mapped parameter `P^3`... 


# Setting non-Gaussian priors

One may also want to impose a specific non-Gaussian prior to some parameters. For example, one may want to set [Beta](https://en.wikipedia.org/wiki/Beta_distribution) or [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) priors on model parameters. Although, superficially, VBA only deals with Gaussian priors, one can use parameter transformations to emulate non-Gaussian priors. This is also exemplified below.

## 

## Emulating sparsity priors

So-called [sparse estimators](https://en.wikipedia.org/wiki/Sparse_approximation) arise in the context of model fitting, when one a priori assumes that only a few (unknown) model parameters deviate from zero. Sparsity constraints can be useful when the estimation problem is under-determined, i.e. when number of model parameters is much higher than the number of data points. In the context of [regularization](https://en.wikipedia.org/wiki/Regularization_(mathematics)) approaches, such constraints are enforced by minimizing [L1](http://mathworld.wolfram.com/L1-Norm.html) or [L0](https://en.wikipedia.org/wiki/Lp_space#When_p_=_0) norms, which yield the so-called [LASSO estimator](https://en.wikipedia.org/wiki/Lasso_(statistics)). 

Note that L1-norm minimization aprpoaches can be seen as a subcase of bayesian parameter estimation under non-Gaussian priors, more precisely: [Laplacian](https://en.wikipedia.org/wiki/Laplace_distribution) priors. In fact, this intuition generalizes to most sparsity constraints (cf. $$Lp$$-norm with $$p<2$$), which have their "sparsity prior" equivalent. It turns out that one can use the following simple parameter transform $$g_s$$ to emulate such sparsity priors without sacrificing the simplicity and robustness of VBA under Gaussian priors:

$$g_s(\vartheta)= (2 s(\vartheta) -1)\vartheta^2$$

where $$s$$ is the standard [sigmoid function](https://en.wikipedia.org/wiki/Sigmoid_function).

This "sparsify mapping" trick is demonstrated in this [technical note](https://arxiv.org/abs/1703.07168), which discloses the properties of such sparse estimators.


# Recovering the re-mapped posterior density

Finally, having performed a VBA analysis with some specified parameter transform in place, one may need to derive summary statistics for the transformed parameters (as opposed to directly using VBA's posterior pdf on native Gaussian parameters). This is the topic of the last section of this page...

## Laplace method

If one is interested in the posterior density over transformed parameters, one can simply use the [Laplace approximation](https://en.wikipedia.org/wiki/Laplace's_method). This is base on a first-order Taylor expansion, which we summarize below.

Let $$f(\vartheta)$$ be the mapping used for setting a given hard constraint on some model parameter $$P=f(\vartheta)$$. Then:
$$ f(\vartheta) = f(E[\vartheta]) + f'(E[\vartheta])\times\left(\vartheta-E[\vartheta]\right) + ... $$
This 




(see `VBA_Laplace.m`).

## Monte-Carlo methods

 Alternatively, one can [sample](https://en.wikipedia.org/wiki/Monte_Carlo_method) (see `VBA_sample.m`) from the Gaussian posterior density over un-transformed parameters, pass the samples through the transform, and then report [summary statistics](https://en.wikipedia.org/wiki/Summary_statistics) over the set of transformed samples (such as mean and variance) and/or full sampling histograms. Such sampling approach can be used, for example, in the aim of recovering [credible intervals](https://en.wikipedia.org/wiki/Credible_interval) over constrained (mapped) parameters.




