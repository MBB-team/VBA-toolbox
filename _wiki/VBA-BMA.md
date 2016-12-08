---
title: "Bayesian model averaging (BMA)"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}



Standard statistical practice ignores model uncertainty. Data analysts typically select a model from some class of models and then proceed as if the selected model had generated the data. This approach ignores the uncertainty in model selection, leading to over-confident inferences. Bayesian model averaging (BMA) provides a coherent mechanism for accounting for this model uncertainty when deriving parameter estimates.

In brief, BMA marginalizes over models to derive posterior densities on model parameters that account for model uncertainty, as follows:
$$p(\theta\mid y) = \sum_{m_i} p(m_i\mid y) p(\theta\mid y,m_i)$$

where $$m_i$$ are the set of candidate models. 
