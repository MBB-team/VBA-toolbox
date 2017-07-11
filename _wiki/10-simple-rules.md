---
title: "Ten simple rules for VBA"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Although this section focuses on good practice with VBA, almost all following "simple rules" actually apply to *any* model-based data analysis. 

# 1: Simulate your model(s)

Most computational models cannot be fully understood from the mathematical equations they entail. So how do you know whether you *know* your model? This is a quote from [Palminteri et al. (2016)](http://www.biorxiv.org/content/early/2016/10/07/079798):

```
The importance of simulating candidate models has been (...) largely overlooked, which entails several drawbacks and leads to invalid conclusions. (...) The analysis of model simulations is often necessary to support the specific claims (...) that most of model-based studies make.
```

We could not have said it better.

Model simulations can (and should) be performed before and after data analysis:

- before the analysis, simulations are useful, e.g., to explore what models can and cannot do, under alternative parameter settings.
- after the analysis, simulations are useful, e.g., to check whether models actually predict specific effects of interest.

> **Tip:** VBA is endowed with specific in-built functions for simulating generative models of any sort. In particular, you can simulate data time series using `simulateNLSS.m` or `simulateNLSS_fb.m` (when on-line feedback is required for simulating the system) and eyeball them using `displaySimulations.m`. Also, in the context of strong non-linearities, full-fledged Monte-Carlo simulations can be derived with `get_MCMC_predictiveDensity.m` and eyeballed using `plotDensity.m`.

Anyway, this is *the* golden rule. Know. Your. Model.


# 2: Start with simple models

The temptation is strong to throw at the data one's most sophisticated model. But this may not be a good idea. First, if you are an unexperienced user of VBA, you may make mistakes (bugs), which may be eventually difficult to correct. Second, your favorite model will benefit from a comparison with simpler models, which will serve as reference points. Third, maybe a simple model will do?


# 3: Perform confusion analyses

VBA offers the possibility of performing principled (Bayesian) model comparisons. Nevertheless, you may not know in advance whether your experimental design enables discriminating between candidate models. Or do you? What about simulating data under alternative models, and then comparing all models given each simulated dataset? This is the essence of *[confusion analyses](https://en.wikipedia.org/wiki/Confusion_matrix)*, which aim at quantifying the expected model selection error rates (under your experimental design). We see this as yet another application of rule #1 :) For an example of how to perform and use confusion analyses, please refer to, e.g., the supplementary materials of [Devaine et al. (2014)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003992).



# 4: Check model identifiability

Interesting models typically include a few unknown parameters, which can be adjusted to observed data in the aim of capturing, e.g., inter-individual differences or treatment effects. But do you know whether the impact of model parameters on predicted data is unambiguous? Can two parameters (or more) predict similar changes in the data? If yes, then the model is not (perfectly) [identifiable](https://en.wikipedia.org/wiki/Identifiability). This is an issue if you are willing to interpret parameter estimates. Thus, in addition to knowing your model, you should know whether it is identifiable, i.e. whether you can recover its parameters from experimental data. This, again, can be assessed using numerical simulations (cf. rule #1).



# 5: Assess the susceptibility of posterior w.r.t. priors

Does posterior inference critically depend upon the priors? If yes, this is a sign that the data is not very informative. In other words, your data do not afford your model. Too bad. Recall that if the data was informative enough, then the posterior would be dominated by the likelihood, i.e. it would be insensitive to the priors. In brief, you may want to perform a [sensitivity analysis](https://en.wikipedia.org/wiki/Sensitivity_analysis), in the aim of evaluating whether you reach the same (posterior) conclusion when using different priors...

> **Tip:** In case posterior inference is highly susceptible to the priors, you can resort to [Bayesian model averaging]({{ site.baseurl }}/wiki/VBA-BMA) or BMA. Recall that different priors effectively induce different generative models. One would then perform model inversion under different priors, and then use BMA to derive posterior estimates that are marginalized over these priors. 



# 6: Optimize your experimental design

Yes, this is a trivial advice. Just in case: did you have a look [here]({{ site.baseurl }}/wiki/Optimizing-the-experimental-design)?


# 7: Derive analytical gradients of evolution/observation functions

VBA is a grown-up: it can work out numerical gradients if you are too lazy to endow evolution and/or observation functions with analytical gradients. But let's face it: by providing these gradients, you could speed up VBA model inversions by orders of magnitude...

> **Tip:** VBA can help you to check whether your derivations are correct by comparing the analytical gradients with the numerical ones. This can be done by setting `options.checkGrads = 1`.



# 8: Update your version of the VBA toolbox regularly

We are making changes on VBA's contents regularly. The only way to keep up is to update your version of the toolbox as regularly as we do. Luckily, [GitHub](https://github.com/MBB-team/VBA-toolbox) can do this automatically for you :)



# 9: Post comments and questions on the forum

Most questions and comments come up multiple times, simply because most people face similar issues. For us, a quick and efficient way of dealing with recurring issues is to keep track of them using the [forum](http://mbb-team.github.io/VBA-toolbox/forum/). Note that many recent VBA developments were motivated by such recurring issues and/or comments...


# 10: Contribute!

You developed a novel computational model for behavioural and/or neurobiological data? You designed an efficient pipeline analysis with existing VBA features? You wrote a cool VBA-compatible piece of code? [Contribute](http://mbb-team.github.io/VBA-toolbox/about/)!
