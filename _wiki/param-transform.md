---
title: "Setting hard constraints through parameter transformations"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


In VBA, priors are Gaussian, i.e. parameters are allowed to vary without bounds (i.e. on the $$\left]-\infty,\infty \right[$$ interval). Setting informative priors allows one to insert "soft constraints" on parameter likely ranges. In particular, decreasing the prior variance eventually insures that parameter posterior estimates will be close to their prior mean. This is called "[shrinkage to the mean](https://en.wikipedia.org/wiki/Shrinkage_estimator)". This sort of constraint is "soft" because there always exist a certain amount of data that eventually allows posterior inference to escape from the influence of the priors. In other terms, as long as the prior variance is not zero, the strength of the prior constraint depends upon the informativeness of the data.

However, one may wan to insert "[hard constraints](https://en.wikipedia.org/wiki/Constrained_optimization)", i.e. some models may include parameters that ought to be bounded within certain ranges. For example, [reinforcement learning](https://en.wikipedia.org/wiki/Reinforcement_learning) (RL) models are parameterized in terms of learning rates, which lie on the $$\left[0,1\right]$$ interval. Thus, to keep the parameters within bounds, one has to [re-parameterize](https://en.wikipedia.org/wiki/Parametrization) the model. Below are two overloaded examples of typical re-parameterizations:

- **range constraint**: one may use the (re-scaled) [sigmoidal](https://en.wikipedia.org/wiki/Sigmoid_function) transform. For example, within a RL evolution function, the learning rate may be defined as `sigm(P)`, where `P` is VBA's evolution parameter (which has a Gaussian prior). Thus, irrespective of the actual value of `P`, the learning rate `sigm(P)` will be constrained to lie on the $$\left[0,1\right]$$ interval.

- **positivity constraint**: one may use the [exponential](https://en.wikipedia.org/wiki/Exponential_function) transform, i.e. `exp(P)` will be positive, irrespective of the actual value of `P`. Another, less "stiff", positivity-enforcing mapping is `exp(1+log(x))`, which behaves linearly away from the origin.

If one is interested in the posterior density over transformed parameters, one can simply use the [Laplace approximation](https://en.wikipedia.org/wiki/Laplace's_method) (see `VBA_Laplace.m`). Alternatively, one can [sample](https://en.wikipedia.org/wiki/Monte_Carlo_method) (see `VBA_sample.m`) from the Gaussian posterior density over un-transformed parameters, pass the samples through the transform, and then report [summary statistics](https://en.wikipedia.org/wiki/Summary_statistics) over the set of transformed samples (such as mean and variance) and/or full sampling histograms. Such sampling approach can be used, for example, in the aim of recovering [credible intervals](https://en.wikipedia.org/wiki/Credible_interval) over constrained (mapped) parameters.


> Note that Bayesian inference is **not invariant through changes in model parameterization**. For example, let us consider the mapping $$P \rightarrow P^3$$. This mapping does not insert any hard constraint, in the sense that the mapped parameter $$P^3$$ is allowed to vary without finite bounds. Nevertheless, Bayesian inference on a model that uses parameter `P` does not yield the same result as inference on the same model, but this time using the mapped parameter `P^3`... 


