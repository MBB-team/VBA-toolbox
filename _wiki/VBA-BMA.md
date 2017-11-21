---
title: "Bayesian model averaging (BMA)"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}



Standard statistical practice ignores model uncertainty. Data analysts typically select a model from some class of models and then proceed as if the selected model had generated the data. This approach ignores the uncertainty in model selection, leading to [over-confident](https://en.wikipedia.org/wiki/Overconfidence_effect) inferences. Bayesian model averaging (BMA) provides a coherent mechanism for accounting for this model uncertainty when deriving parameter estimates.

In brief, BMA marginalizes over models to derive posterior densities on model parameters that account for model uncertainty, as follows:

$$p(\theta\mid y) = \sum_{m_i} p(m_i\mid y) p(\theta\mid y,m_i)$$

where $$m_i$$ are the set of candidate models, $$p(m_i\mid y)$$ is the posterior probability over model $$m_i$$, and $$p(\theta\mid y,m_i)$$ is the posterior density on model parameters conditional on model $$m_i$$. The latter posterior density is a decent proxy for one's information on parameters $$\theta$$ only if $$p(m_i\mid y) \approx 1$$. Otherwise, uncertainty regarding the correct model will automatically translate in uncertainty regarding model parameters...


In VBA, BMA can be performed by first inverting the candidate models, and then calling the BMA routine:

`[p_BMA] = VBA_BMA(p0,F0)`.

whose inputs and outputs are described below:

- `p0`: a Kx1 cell-array of VBA posterior structures, which are conditional onto specific [generative models](https://en.wikipedia.org/wiki/Generative_model) (where K is the number of candidate models)
- `F0`: a Kx1 vector of log-[model evidences](https://en.wikipedia.org/wiki/Marginal_likelihood) (for each candidate model)
- `p_BMA`: the resulting posterior structure, with the first two [moments](https://en.wikipedia.org/wiki/Moment_(mathematics)) of the marginal probability density functions.

Note that `VBA_BMA.m` provides a Gaussian approximation of the marginal posterior over parameters.
