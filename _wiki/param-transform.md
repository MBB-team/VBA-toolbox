---
title: "Setting hard constraints through parameter transformations"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


In VBA, priors are Gaussian, i.e. parameters are allowed to vary without bounds (i.e. on the ]-Inf,Inf[ interval). Setting informative priors allows one to insert "soft constraints" on parameter likely ranges. In particular, decreasing the prior variance eventually insures that parameter posterior estimates will be close to their prior mean. This is called "shrinkage to the mean". This sort of constraint is "soft" because there always exist a certain amount of data that eventually allows posterior inference to escape from the influence of the priors. In other terms, as long as the prior variance is not zero, the strength of the prior constraint depends upon the informativeness of the data.

However, one may wan to insert "hard constraints", i.e. some models may include parameters that ought to be bounded within certain ranges. For example, reinforcement learning models are parameterized in terms of learning rates, which lie on the [0,1] interval. Thus, to keep the parameters within bounds, one has to re-parameterize the model. Below are two overloaded examples of typical re-parameterizations:

- **range constraint**: one may use the (re-scaled) sigmoidal transform. For example, within a RL evolution function, the learning rate may be defined as `sigm(P)`, where `P` is VBA's evolution parameter (which has a Gaussian prior). Thus, irrespective of the actual value of `P`, the learing rate `sigm(P)` will be constrained to lie on the [0,1] interval.

- **positivity constraint**: one may use the exponential transform, i.e. `exp(P)` will be positive, irrespective of the actual value of `P`. 

If one is interested in the posterior density over transformed parameters, one can simply use the Laplace approximation (see `VBA_Laplace.m`). Alternatively, one can sample (see `VBA_sample.m`) from the Gaussian posterior density over un-transformed parameters, pass the samples through the transform, and then report summary statistics over the set of transformed samples (such as mean and variance) and/or full sampling histograms. 

Note that Bayesian inference is not invariant through changes in model parameterization. For example, let us consider the mapping $$P \rightarrow P^3$$. This mapping does not insert any hard constraint, in the sense that the mapped parameter $$P^3$$ is allowed to vary without finite bounds. Inference on a model that uses parameter `P` does not yield the same result as inference on the same model, but this time using the mapped parameter `P^3`... 


