---
title: "Setting hard constraints through parameter transformations"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


In VBA, priors are Gaussian, i.e. parameters are allowed to vary without bounds (i.e. on the $$\left]-\infty,\infty \right[$$ interval). Setting informative priors allows one to insert "soft constraints" on parameter likely ranges. In particular, decreasing the prior variance eventually insures that parameter posterior estimates will be close to their prior mean. This is called "[shrinkage to the mean](https://en.wikipedia.org/wiki/Shrinkage_estimator)". This sort of constraint is "soft" because there always exist a certain amount of data that eventually allows posterior inference to escape from the influence of the priors. In other terms, as long as the prior variance is not zero, the strength of the prior constraint depends upon the informativeness of the data.

However, one may want to insert "[hard constraints](https://en.wikipedia.org/wiki/Constrained_optimization)", i.e. some models may include parameters that ought to be bounded within certain ranges. For example, [reinforcement learning](https://en.wikipedia.org/wiki/Reinforcement_learning) (RL) models are parameterized in terms of learning rates, which lie on the $$\left[0,1\right]$$ interval. Thus, to keep the parameters within bounds, one has to [re-parameterize](https://en.wikipedia.org/wiki/Parametrization) the model. Below are two overloaded examples of typical re-parameterizations for setting such hard constraints.



# Examples of re-parameterizations for setting hard constraints

## Range constraint

One may use the (re-scaled) [sigmoidal](https://en.wikipedia.org/wiki/Sigmoid_function) transform $$s(x) = \frac{1}{1+e^{-x}}$$. For example, within a RL evolution function, the learning rate $$z$$ may be defined as $$z=s(x)$$, where $$x$$ is VBA's evolution parameter (which has a Gaussian prior). Thus, irrespective of the actual value of $$x$$, the learning rate $$z$$ will be constrained to lie on the $$\left[0,1\right]$$ interval.

## Positivity constraint

One may use the [exponential](https://en.wikipedia.org/wiki/Exponential_function) transform, i.e. $$e^x$$ will be positive, irrespective of the actual value of $$x$$. Another, less "stiff", positivity-enforcing mapping is $$e^{1+log(x)}$$, which behaves like the identity mapping when $$x$$ is sufficiently far from the origin.


> Note that Bayesian inference is **not invariant through changes in model parameterization**. For example, let us consider the mapping $$x \rightarrow x^3$$. This mapping does not insert any hard constraint, in the sense that the mapped parameter $$x^3$$ is allowed to vary without finite bounds. Nevertheless, Bayesian inference on a model that uses parameter $$x$$ does not yield the same result as inference on the same model, but this time using the mapped parameter $$x^3$$. This is because passing through the $$x \rightarrow x^3$$ mapping can be seen as setting a non-Gaussian prior on $$x^3$$...


# Setting non-Gaussian priors

One may also want to impose a specific non-Gaussian prior to some parameters. For example, one may want to set [Beta](https://en.wikipedia.org/wiki/Beta_distribution) or [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) priors on model parameters. Although, superficially, VBA only deals with Gaussian priors, one can use parameter transformations to emulate non-Gaussian priors. This is exemplified below.

## Setting a specific prior density

Let $$z$$ be some model parameter, which one may want to set some specific prior density $$p_z(z)$$ on (e.g. Beta, Gamma, etc...). This can be done by mapping VBA's native parameter $$x$$ through the following transformation $$h(x)$$ (in the corresponding evolution or observation function):

$$ h(x)=P_z^{-1}\left(\Phi\left(x\right)\right) $$

where $$\Phi(x)$$ is the cumulative distribution function of VBA's Gaussian priors on $$x$$ and $$P_z^{-1}\left(z\right)$$ is the inverse cumulative distribution function of $$z$$. The transformation $$h(x)$$ is the composition of two mappings:
- the [probability integral transform](https://en.wikipedia.org/wiki/Probability_integral_transform), which produces a uniform distribution over [0,1] from any (here: gaussian) density,
- and the [inverse transform sampling](https://en.wikipedia.org/wiki/Inverse_transform_sampling), which produces any distribution (here: $$p_z(z)$$) from the uniform distribution over [0,1].

The ensuing transformed parameter $$z=h(x)$$ follows the distribution $$p_z(z)$$, irrespective of the way VBA's native prior is specified, as long as the tansformation $$h(x)$$ uses the corresponding normal cumulative distribution function (with identical mean and variance). 

> NB: this section was inspired from a comment made by **Emma** on VBA's forum: thank you!


## Emulating sparsity priors

So-called [sparse estimators](https://en.wikipedia.org/wiki/Sparse_approximation) arise in the context of model fitting, when one a priori assumes that only a few (unknown) model parameters deviate from zero. Sparsity constraints can be useful when the estimation problem is under-determined, i.e. when the number of model parameters is much higher than the number of data points. In the context of [regularization](https://en.wikipedia.org/wiki/Regularization_(mathematics)) approaches, such constraints are enforced by minimizing [L1](http://mathworld.wolfram.com/L1-Norm.html) or [L0](https://en.wikipedia.org/wiki/Lp_space#When_p_=_0) norms, which yield the so-called [LASSO estimator](https://en.wikipedia.org/wiki/Lasso_(statistics)). 

Note that L1-norm minimization approaches can be seen as a subcase of bayesian parameter estimation under non-Gaussian priors, more precisely: [Laplacian](https://en.wikipedia.org/wiki/Laplace_distribution) priors. In fact, this intuition generalizes to most sparsity constraints (cf. $$Lp$$-norm with $$p<2$$), which have their "sparsity prior" equivalent. For example, emulating L1-like sparsity priors can be done by mapping VBA's native parameters $$x$$ through the following simple transform: 

$$ g_s(x)= \left(2 s(x) -1\right)x^2 $$

where $$s(x)$$ is the standard sigmoid function.

> This trick is justified and demonstrated in this [technical note](https://arxiv.org/abs/1703.07168), which discloses the properties of the ensuing (sparse) VBA estimators. Note that, to ensure that VBA's estimation converges, one has to ensure that the prior mean of native VBA parameters $$x$$ is not exactly equal to 0. This is because otherwise, the gradient of $$g_s(x)$$ w.r.t. model parameters is null, which prevents VBA from updating the prior.


# Recovering the re-mapped posterior density

Finally, having performed a VBA analysis where native parameters $$x$$ have been mapped through some specific transform $$g(x)$$, one may need to derive summary statistics for the transformed parameters $$z=g(x)$$ (as opposed to directly using VBA's posterior pdf on native Gaussian parameters). This is the topic of the last section of this page.

## A few specific cases

We begin by highlighting simple analytic expressions for the first- and second-order moments of gaussian variables passed though the exponential (cf. positivity constraint) and sigmoidal (cf. range constraint) mappings.

### The exponential mapping

This is useful for imposing positivity constraints on model parameters (see above). Let $$g(x)=e^x$$ and $$x$$ be normally distributed, i.e.: $$p(x) = N\left(\mu,\sigma\right)$$. Then:

$$ E[g(x)] = e^{\mu + \frac{\sigma}{2}} $$,

and

$$ V[g(x)] = e^{2\mu + \sigma} \left(e^{\sigma}-1\right) $$,

Recovering the two first two moments of exponentially-mapped parameters $$z=e^x$$ can thus be done analytically, given the first two moments of the Gaussian posterior distribution on $$x$$.

### The sigmoid mapping

This is useful for imposing "range" constraints on model parameters (see above). Let $$s(x)=\frac{1}{1+e^{-x}}$$ and $$x$$ be normally distributed, i.e.: $$p(x) = N\left(\mu,\sigma\right)$$. Then:

$$ E[s(x)] \approx s\left(\frac{\mu}{\sqrt{1+a\sigma}}\right) $$,

and

$$ V[s(x)] \approx s\left(\frac{\mu}{\sqrt{1+a\sigma}}\right)\left(1-s\left(\frac{\mu}{\sqrt{1+a\sigma}}\right)\right)\left(1-\frac{1}{\sqrt{1+a\sigma}}\right) $$,

where $$a=\frac{3}{\pi^2} \approx 0.3$$.

> These approximations, and other related ones (e.g., softmax with more than 2 variables, log-sigmoid transforms, etc...) are described in [this technical note](https://arxiv.org/abs/1703.00091). Suffices to say here that they yield less than 2% relative error.



## Laplace's method

If no analytical approximation exists for the moments of a Gaussian variable $$x$$ passed through the mapping $$g(x)$$, one can use the so-called [Laplace approximation](https://en.wikipedia.org/wiki/Laplace's_method), which we summarize below.

Let $$g(x)$$ be the mapping used for setting a given hard constraint on some model parameter $$z=g(x)$$. Then a (truncated) first-order [Taylor expansion](https://en.wikipedia.org/wiki/Taylor_series) in the vicinity of $$E[x]$$ yields:

$$ g(x) = g(E[x]) + g'(E[x])\times\left(x-E[x]\right) + ... $$ 

This first-order Taylor expansion can be used to derive the first- and second-order moments of $$z$$, given first- and second-order moments of $$x$$:

$$ E[z] \approx g(E[x]) $$

and

$$ V[z] \approx V[x]\times g'(E[x])^2 $$
 
The function `VBA_getLaplace` can be used to derive the above Laplace approximation, as follows:

```matlab
dim = struct('n',0,'n_theta',0,'n_phi',1);
opt.priors.muPhi = posterior.muTheta(ind);
opt.priors.SigmaPhi = posterior.SigmaTheta(ind,ind);
[Ez,Vz] = VBA_getLaplace([],[],@myMapping,dim,opt,0);
```

where `ind` is the index of the evolution parameter that went through the transform, `@myMapping` implements the parameter tansformation $$g$$ (but with the usual i/o format of VBA observation functions), `posterior` has been obtained using VBA, and `Ez` and `Vz` are the Laplace approximations to the first- and second-order moments of $$z$$, respectively...


## Monte-Carlo's method

Alternatively, one can [sample](https://en.wikipedia.org/wiki/Monte_Carlo_method) from the Gaussian posterior density over un-transformed parameters, pass the samples through the transform, and then report [summary statistics](https://en.wikipedia.org/wiki/Summary_statistics) over the set of transformed samples (such as mean and variance). This is the essence of the [Monte-Carlo method](https://en.wikipedia.org/wiki/Monte_Carlo_method), which are widely used in statistics and numerical analysis. 

The following lines of code reproduce Monte-Carlo's method for the same example as the section above:

```matlab
suffStat.mu = posterior.muTheta(ind);
suffStat.Sigma = posterior.SigmaTheta(ind,ind);
N = 1e6; % number of Monte-Carlo samples
X = VBA_sample('gaussian',suffStat,N,0); % sample from gaussian posterior 
gX = feval(@myMapping,[],X,[],[]); % cf. VBA i/o structure
Ez = mean(gX);
Vz = var(gX);
```

One possible advantage of Monte-Carlo's method is that one can report more than just the mean and variance of the mapped parameters. In particular, one can derive [credible intervals](https://en.wikipedia.org/wiki/Credible_interval) over constrained (mapped) parameters, and/or full sampling histograms.

> The accuracy of Monte-Carlo's method depends upon the numer of samples that are used to derive the moments of the distribution (the rule of thumb is typically $$10^5$$ samples at least for each degree of freedom). This eventually induces a very high computational cost when the dimensionality of the parameter space increases...

