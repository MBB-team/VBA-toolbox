---
title: "Ten simple rules for good practice with VBA"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

# 1: Simulate observable data under your model(s)

This is a quote from [Palminteri et al. (2016)](http://www.biorxiv.org/content/early/2016/10/07/079798):

```
Cognitive neuroscience, especially in the fields of learning and decision-making, is witnessing the blossoming of computational model-based analyses. [...] Candidate models should be compared by trading off their ability to predict the data as a function of their complexity. However, the importance of simulating candidate models has been so far largely overlooked, which entails several drawbacks and leads to invalid conclusions. Here we argue that the analysis of model simulations is often necessary to support the specific claims about behavioral function that most of model-based studies make.
```

We could not have said it better.

Model simulations should be performed before and after data analysis:

- before the analysis, simulations are useful to explore what models can and cannot do, given alternative parameter settings.
- after the analysis, simulations are useful to predict yet unobserved data (maybe under different experimental conditions).

This **the ** golden rule. Know your model(s).


# 2: Start with simple models

Take our word: the temptation is strong to throw at the data one's most sophisticated model. But this may not be a good idea. First, if you are an unexperienced user of VBA, you may make mistakes (bugs), which may be eventually difficult to correct. Second, your favorite model will benefit from a comparison with simpler models, which will serve as reference points. Third, maybe, who knows, a simple model will do?


# 3: Perform confusion analyses

VBA offers the possibility of performing principled (Bayesian) model comparisons. Nevertheless, you may not know in advance whether your data yields enough statistical power for discriminating between alternative models. Or do you? What about simulating data under different models, and then compare all models given each simulated dataset? This is the essence of *confusion analyses*, which aim at quantifying the expected model selection error rates (under your expeirmental design). For exmaple of how to do this, please see the supplmentary materials of, e.g. [Devaine et al. (2014)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003992).

> This is yet another application of the rule #1...


# 4: Check model identifiability





# 4: Derive analytical gradients of evolution/observation functions

# 5: Assess the susceptibility of posterior w.r.t. priors


# 7: Optimize your experimental design

# 8: Update your version of the VBA toolbox regularly

# 9: Post comments and questions on the forum

# 10: Contribute!
