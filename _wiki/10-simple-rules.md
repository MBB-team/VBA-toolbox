---
title: "Ten simple rules for VBA"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Although this section focuses on good practice with VBA, almost all following "simple rules" actually apply to *any* model-based data analysis. 

# 1: Simulate your model(s)

This is a quote from [Palminteri et al. (2016)](http://www.biorxiv.org/content/early/2016/10/07/079798):

```
Candidate models should be compared by trading off their ability to predict the data as a function of their complexity. However, the importance of simulating candidate models has been so far largely overlooked, which entails several drawbacks and leads to invalid conclusions. Here we argue that the analysis of model simulations is often necessary to support the specific claims (...) that most of model-based studies make.

```

We could not have said it better.

Model simulations should be performed before and after data analysis:

- before the analysis, simulations are useful to explore what models can and cannot do, given alternative parameter settings.
- after the analysis, simulations are useful to predict yet unobserved data (maybe under different experimental conditions), for experimental validation purposes.

This **the** golden rule. Know your model(s).


# 2: Start with simple models

Take our word: the temptation is strong to throw at the data one's most sophisticated model. But this may not be a good idea. First, if you are an unexperienced user of VBA, you may make mistakes (bugs), which may be eventually difficult to correct. Second, your favorite model will benefit from a comparison with simpler models, which will serve as reference points. Third, maybe, who knows, a simple model will do?


# 3: Perform confusion analyses

VBA offers the possibility of performing principled (Bayesian) model comparisons. Nevertheless, you may not know in advance whether your data yields enough statistical power for discriminating between alternative models. Or do you? What about simulating data under different models, and then compare all models given each simulated dataset? This is the essence of *confusion analyses*, which aim at quantifying the expected model selection error rates (under your experimental design). We see this as yet another application of rule #1 :) For an example of how to do this, please see the supplementary materials of, e.g. [Devaine et al. (2014)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003992).



# 4: Check model identifiability

Interesting models typically include a few unknown parameters, which can be adjusted to observed data in the aim of capturing, e.g., inter-individual differences or treatment effects. But do you know whether the impact of model parameters on predicted data is unambiguous? Can two parameters (or more) predict similar changes in the data? If yes, then the model is not (perfectly) identifiable. This may be an issue, if you are willing to interpret these parameters. Thus, in addition to knowing your model, you should know whether it is identifiable, i.e. whether you can recover its parameters of interest from experimental data. This, again, can be assessed using numerical simulations (rule #1).



# 5: Assess the susceptibility of posterior w.r.t. priors

A related issue is whether posterior inference critically depends upon the priors. This, in itself, is a sign that the data is not informative enough. Too bad. Recall that if the data was highly informative, then the posterior would be dominated by the likelihood, i.e. it would be insensitive to the priors. Thus, you may want to know whether you reach the same (posterior) conclusion when using different priors...



# 6: Optimize your experimental design

Yes, this is a trivial advice. No comment (did you have a look [here]({{ site.baseurl }}/wiki/Optimizing-the-experimental-design)?).


# 7: Derive analytical gradients of evolution/observation functions

VBA is a grown-up: it can work out numerical gradients if you are too lazy to endow evolution and/or observation functions with analytical gradients. But let's face it: by providing these gradients, you could speed up VBA model inversions by orders of magnitude. This is perhaps the simplest and most efficient way of improving your code...

> **Tip:** VBA can help you to check whether your derivations are correct by comparing the analytical gradients with the numerical ones. This can be done by setting `options.checkGrads = 1`.



# 8: Update your version of the VBA toolbox regularly

We are making changes on VBA's content regularly. The only way to keep up is to update your version of the toolbox as regularly as we do. Luckily, [GitHub](https://github.com/MBB-team/VBA-toolbox) can do this automatically for you :)



# 9: Post comments and questions on the forum

Most questions and comments come up multiple times, simply because most people face the same sorts of issues. For us, a quick and efficient way of dealing with recurring issues is to keep track of them using the [forum](http://mbb-team.github.io/VBA-toolbox/forum/)...


# 10: Contribute!

You developed an exciting computational model for behavioural and/or neurobiological data? You designed an efficient pipeline analysis with existing VBA features? You wrote a cool VBA-compatible piece of code? [Contribute](http://mbb-team.github.io/VBA-toolbox/about/)!
