---
title: "Bayesian modeling: introduction"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

# Derivation of the likelihood function

One usually starts with a quantitative assumption or model of how observations y are generated. Without loss of generality, this model possesses unknown parameters $$\vartheta$$, which are mapped through an observation function g:

$$y= g(\vartheta)+\epsilon$$

where ε are model residuals or measurement noise. If the (physical) processes underlying ε were known, they would be included in the deterministic part of the model, i.e.: $$\varepsilon \in g(\vartheta)$$. Typically, we thus have to place statistical priors on $$\varepsilon$$, which eventually convey our lack of knowledge, as in “the noise is small”. This can be formalized as a probabilistic statement, such as: “the probability of observing big noise is small”. Under the central limit theorem, such prior would be equivalent to assuming the noise follows a normal distribution:

$$p(\varepsilon\mid m)\propto exp\left(-\frac{1}{2\sigma^2}\varepsilon^2\right) \implies p(\lvert\epsilon\rvert>1.96\sigma\mid m) \approx 0.05$$

where $$\sigma$$ is the noise’ standard deviation (it determines how big is “big”) and m is the so-called generative model. Equations 1 and 2 are compiled to derive a likelihood function $$p(y\mid\vartheta,m)$$, which specifies how likely it is to observe any particular set of observations y, given the unknown parameters $$\vartheta$$ of the model m :

$$p(y\mid\vartheta,m) = exp\left(-\frac{1}{2\sigma^2}(y-g(\vartheta))^2\right)$$

The intuition underlying the above derivation of the likelihood function can be generalized to any generative model $$m$$, whose parameters $$\vartheta$$ simply control the statistical moments of the distribution $$p(y\mid\vartheta,m)$$. The key point here is that the likelihood function always derives from priors about observation mappings and measurement noise.

# Bayes' rule

The likelihood function is the statistical construct that is common to both frequentist (classical) and bayesian inference approaches. However, bayesian approaches also require the definition of a prior distribution $$p(\vartheta\mid m)$$ on model parameters $$\vartheta$$, which reflects knowledge about their likely range of values, before having observed the data y. As for priors about measurement noise, such priors can be (i) principled (e.g. certain parameters cannot have negative values), (ii) conservative (e.g. “shrinkage” priors that express the assumption that coupling parameters are small), or (iii) empirical (based on previous, independent measurements).
Combining the priors and the likelihood function allows one, via Bayes' Theorem, to derive both the marginal likelihood of the model (the so-called model evidence):

$$p(y\mid m)=\int p(y\mid \vartheta,m)p(\vartheta\mid m)d\vartheta$$

and the posterior probability density function $$p(\vartheta\mid,m)$$ over model parameters $$\vartheta$$:

$$p(\vartheta\mid,m)=\frac{p(y\mid\vartheta,m)p(\vartheta\mid m)}{p(y\mid m)}$$

This is called “model inversion” or “solving the inverse problem”. The posterior density $$p(\vartheta\mid y, m)$$  quantifies how likely is any value of $$\vartheta$$, given the observed data y and the generative model $$m$$. It is used for inferring on “interesting” model parameters, by marginalizing over “nuisance” parameters. The model evidence $$p(y\mid m)$$  quantifies how likely is the observed data y under the generative model $$m$$. Another perspective on this is that $$-\log p(y\mid m)$$ measures statistical surprise, i.e. how unpredictable was the observed data $$y$$ under the model $$m$$. The model evidence accounts for model complexity, and thus penalizes models, whose predictions do not generalize easily (this is referred to as “Occam’s razor”). Under flat priors over models, it is used for model selection (by comparison with other models that differ in terms of either their likelihood or their prior density).


# Statistical tests and Bayesian model comparison

Moments of the posterior density $$p(\vartheta\mid y,m)$$ can be used to define parameter estimates (e.g., the posterior mean). Typically, as the quantity of available data increases, Bayesian parameter estimates effectively converge to frequentist (e.g. maximum likelihood) estimators. This is because the weight of the prior on any moment of the posterior distribution becomes negligible.

However, this (asymptotic) equivalence does not hold for model comparison. This is important, because model comparison has many application within a Bayesian framework. For example, when testing whether a parameter is zero, one effectively compares two hypotheses: the 'null', in which the parameter is fixed to zero, against the 'alternative', in which the parameter is allowed to vary.

According to the Neyman-Pearson lemma, the most powerful test to compare such two hypotheses or models is the likelihood-ratio test, i.e.:

$$\frac{p(y\mid m_1)}{p(y\mid m_2)} > K$$

where $$K$$ is set to satisfy a controlled statistical risk.
This motivates the use of model evidences to perform statistical testing (e.g. testing the null) within a Bayesian framework. In fact, the quantity above is known as the 'Bayes factor', and is used whenever one wants to select between two models. More generally, the comparison of more than two models can be based upon their model evidence, which quantifies the plausibility of the data under any model.

The critical thing to note is that the model evidence is not a simple measure of model fit: there is an inherent penalization for model complexity. Let us note that the model evidence is maximal whenever the information in the data about model parameters is the same as that of the prior density. In other words, a good model is a model that is confirmed by the data. A simple model has tight priors: to the limit, the simplest of all models has no unknown parameters (infinite precision priors). This means that there is a complexity cost to changing one's belief on parameters to explain data.

In brief, the model evidence is nothing else than the data likelihood, given any model. Bayes' rule can then be used to perform inference on models, by deriving the posterior distribution over models, i.e.:

$$p(m\mid y)=\frac{p(y\mid m)p(m)}{p(y)}$$

where $$p(y)$$ is the probability of data given all possible models

$$p(y)=\sum _m p(y\mid m)p(m)$$

The term $$p(m)$$ is the prior probability on model m. Typically, noninformative priors are used and the equation above is driven solely by the Bayes' factor.

As for any statistical test, a threshold has to be set for deciding whether a model is "better" than another one. This threshold can be chosen similarly to classical statistics, i.e. on the basis of some acceptable statistical risk. It turns out that the probability of making a model selection error is 1 minus the posterior probability of the selected model. If this probability has to be controlled at e.g., 0.05, then one "selects" a model only if its posterior probability exceeds 0.95. When comparing two models with each other, this corresponds to a threshold of $$K=20$$ on the Bayes' factor.

This logic applies for any set of models, given any data. The only constraint is that model evidences have to be evaluated on the same data  set.

Practically speaking, the Bayes' factor induces three types of statistical decisions:

- $$K>20$$:      model $$m_1$$ is selected
- $$0.05<K<20$$: no model is selected
- $$K<0.05$$:    model $$m_2$$ is selected
