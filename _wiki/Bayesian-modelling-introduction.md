---
title: "Bayesian inference: introduction"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

# Derivation of the likelihood function

One usually starts with a quantitative assumption or model of how observations y are generated. Without loss of generality, this model possesses unknown parameters $$\vartheta$$, which are mapped through an observation function $$g$$:

$$y= g(\vartheta)+\epsilon$$

where $$\epsilon$$ are model residuals or measurement noise. If the (physical) processes underlying $$\epsilon$$ were known, they would be included in the deterministic part of the model. Typically, we thus have to place statistical priors on $$\epsilon$$, which eventually convey our lack of knowledge, as in “the noise is small”. This can be formalized as a probabilistic statement, such as: “the probability of observing big noise is small”. Under the central limit theorem, such prior would be equivalent to assuming that the noise follows a normal distribution:

$$p(\epsilon\mid m)\propto exp\left(-\frac{1}{2\sigma^2}\varepsilon^2\right) \implies p(\lvert\epsilon\rvert>1.96\sigma\mid m) \approx 0.05$$

where $$\sigma$$ is the noise’ standard deviation (it determines how big is “big”) and $$m$$ is the so-called generative model. Combining the two equations above yields the likelihood function $$p(y\mid\vartheta,m)$$, which specifies how likely it is to observe any particular set of observations $$y$$, given the unknown parameters $$\vartheta$$ of the model $$m$$ :

$$p(y\mid\vartheta,m) = exp\left(-\frac{1}{2\sigma^2}(y-g(\vartheta))^2\right)$$

This derivation of the likelihood function can be generalized to any generative model $$m$$, whose parameters $$\vartheta$$ simply control the statistical moments of the distribution $$p(y\mid\vartheta,m)$$. The key point here is that the likelihood function always derives from prior assumptionss about observation mappings $$g(\vartheta)$$ and measurement noise $$\epsilon$$.

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
This motivates the use of model evidences to perform statistical testing (e.g. testing the null) within a Bayesian framework. In fact, the quantity above is known as the **Bayes' factor**, and is used whenever one wants to select between two models. Practically speaking, the Bayes' factor induces three types of statistical decisions:

- $$K>20$$:      model $$m_1$$ is selected
- $$0.05<K<20$$: no model is selected
- $$K<0.05$$:    model $$m_2$$ is selected


The critical thing to note is that the model evidence $$p\big( y\mid m\big)$$ is not a simple measure of model fit: there is an inherent penalization for model complexity. This penalization is intimately related to the priors. In brief, a simple model has tight priors: at the limit, the simplest model has no unknown parameters (infinite prior precision). More complex models are equipped with vague priors, which will be updated to a larger extent once the data has been observed. However, this flexibility has a cost: that of confusing noise $$\epsilon$$ with variations in the data that are induced by $$g(\vartheta)$$. This is called "over-fitting" the data, and results in greater error when extrapolating model predictions.

Model evidence is essentially a trade-off between fit accuracy and model complexity. It can be used to compare more than two models, simply because it quantifies the plausibility of the data under any model. In fact, Bayesian model comparison proceeds with the exact same logic than when performing inference on parameters. Here again, one relies upon Bayes' rule to derive the posterior distribution over models, i.e.:

$$p(m\mid y)=\frac{p(y\mid m)p(m)}{p(y)}$$

where $$p(y)$$ is the probability of data given all possible models:

$$p(y)=\sum _m p(y\mid m)p(m)$$

The term $$p(m)$$ is the prior probability on model $$m$$. Typically, non-informative priors are used and the equation above is driven solely by model evidences.

As for any statistical test, a threshold has to be set for deciding whether a model is "better" than another one. This threshold can be chosen similarly to classical statistics, i.e. on the basis of some acceptable statistical risk. It turns out that the probability of making a model selection error is 1 minus the posterior probability of the selected model. If this probability has to be controlled at e.g., 0.05, then one "selects" a model only if its posterior probability exceeds 0.95. When comparing two models with each other, this corresponds to a threshold of $$K=20$$ on the Bayes' factor.

This reasoning applies for any set of models, given any data. The only constraint is that model evidences have to be evaluated on the same data  set.


# Classical versus Bayesian hypothesis testing

Let us first summarize how classical (frequentist) testing proceeds. One starts with defining *the null*, e.g., some parameter of interest has no effect ($$H_0: \theta = 0$$). One then construct a test statistic $$t$$ (e.g., Student t-test), for which the distribution under the null $$p\left(t \mid H_0 \right)$$ is known. One then evaluates the test statistic given one's data ($$t^\ast$$) and compare it to the distribution $$p\left(t \mid H_0 \right)$$. Typically, one wants to control the [false positive rate](https://en.wikipedia.org/wiki/False_positive_rate), i.e. the probability of rejecting the null while it is true. This is why the classical decision rule is something like this: if $$P\left(t>t^\ast \mid H_0 \right) < \alpha$$, then reject the null (implicitly: accept the alternative, i.e.: $$H_1: \theta \neq 0$$). This is exmplified in the figure below:

![]({{ site.baseurl }}/images/wiki/pval.jpg)

The distribution under the null $$p\left(t \mid H_0 \right)$$ is shown in red, along with the probability of finding a more extreme value under the null $$P\left(t>t^\ast \mid H_0 \right)$$ (grey area under the cruve). Under this decision rule, we expect to wrongly accept the null in $$100 \times \alpha \%$$ cases.

But what happens when $$P\left(t>t^\ast \mid H_0 \right) > \alpha$$? Is this evidence in favour of $$H_0$$? Of course not. At this point, one has to resort to the notion of [statistical power](https://en.wikipedia.org/wiki/Statistical_power), i.e. the test's sensitivity. Formally speaking, sensitivity is one minus the [false negative rate](https://en.wikipedia.org/wiki/False_positives_and_false_negatives), i.e. the probability of not detecting an effect when there is. In principle, this probability can be derived from the alternative hypothesis $$H_1: \theta \neq 0$$. For example, if it was the case that, for $$\alpha=0.05$$, power was above 95%, then one could formally *accept the null*, because $$P\left(t < t^\ast \mid H_1 \right) < 0.05 $$. This logic, however, is more appropriately pursued using Bayesian hypothesis testing.

Recall that Bayesian hypothesis testing is a special case of Bayesian model comparison, whereby one would compare $$H_0$$ with $$H_1$$. This requires one to consider the alternative hypothesis explicitly. But what novel perspective does it bring? To begin with, although a given data sample may be unlikely under $$H_0$$, it may be even more unlikey under $$H_1$$. In fact, Bayesian model comparison is more concerned with total error rates than with false positive rates alone. Let us extend our exemple above:

![]({{ site.baseurl }}/images/wiki/bmc.jpg)

The distribution of possible data samples $$y$$ under the null ($$p\left(y \mid H_0 \right)$$) and under the alternative ($$p\left(y \mid H_1 \right)$$) are shown in red and in green, respectively. Under uninformative priors on models, Bayesian model comparison would favour $$H_0$$ whenever $$p\left(y \mid H_0 \right) > p\left(y \mid H_1 \right)$$, and $$H_1$$ otherwise. The ensuing Bayesian decision criterion is depicted using the black solid line, which minimizes the total error rate (grey area under the cruve). In this example, Bayesian model comparison would yield more false positives than classical hypothesis testing (cf. dotted black line). However, one can see that this is more than compensated by the higher sensitivity of Bayesian approaches, which yield much smaller false negative rate than classical hypothesis testing.

Note that Bayesian model comparison is not necessarily more liberal than classical hypothesis testing:

![]({{ site.baseurl }}/images/wiki/bmc2.jpg)

In this example, the alternative hypothesis $$H_1$$ is further away from $$H_0$$, and Bayesian model comparison would yield less false positives than clssical hypothesis testing. This example is in fact typical of situations where, although a given data sample may be unlikely under the null (according to classical hypothesis testing), it may still be even more unlikely under the alternative...







