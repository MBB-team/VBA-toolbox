---
title: "Bayesian learning"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Bayesian decision theory (BDT) is a probabilistic framework that is concerned with how decisions are made or should be made, in ambiguous or uncertain situations. The “normative” approach to BDT focuses on identifying the “optimal” or “rational” decision in a given context. This has been applied extensively to methods supporting decisions (e.g., statistical testing). In contradistinction, “descriptive” BDT attempts to describe what people actually do. In this context, optimal decisions are considered as quantitative predictions that can be tested against observed behaviour.

BDT relies on two processes: belief updating and decision making, which are related to the key elements of BDT; namely, prior distributions and utility functions, respectively. In the context of perceptual categorisation, prior beliefs are motivated by the inherent ambiguity of sensory information, which leads to uncertainty about the underlying causes of sensory signals. Priors effectively resolve this ambiguity and are thought to be the basis of most sensory illusions and multistable perceptual effects. In addition, BDT is bound to a perspective on preferences, namely “utility theory”, which was explored in length in the context of economic decisions. In this context, utility functions can be regarded as a surrogate for a task goal or, equivalently, a scoring of the subject’s preferences.

Bayesian models can be used either to capture peoples' hidden priors beliefs (which shape the way they learn) and/or preferences (which map beliefs onto actions).

Many Bayesian models are available in the VBA toolbox.

The script `demo_AVL_recog.m` demonstrates a simple Bayesian learning model, which predicts trial-by-trial variations in reaction times during an associative learning task. The constituents of this model are:

- a set of hierarchically four organized states (3 per available action). The first level captures potential perceptual uncertainty in sensory feedbacks. The second level contains the moments (mean and variance) of peoples' posterior belief about the cue-outcome contingency.
- two evolution parameters that control the perceptual uncertainty and the prior variance of cue-action contingency, respectively.
- two observation parameters that control the speed-accuracy trade-off during perceptual categorization.

We refer the interested reader to [Daunizeau et al. 2010](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0015555).

The script `demo_volatileVB.m` demonstrates a (hierarchical) Bayesian learning model in the context of operant learning (two-armed bandit task). The constituents of this model are:

- a set of hierarchically organized states (6 per available action). The first two levels are similar to the above model. The third level contains the moments (mean and variance) of peoples' posterior belief about the the action-outcome contingency's volatility. Critically, the effective learning rate of such an agent follows his belief about the environmental volatility.
- three evolution parameters that control the dynamical changes of the agent's effective learning rate (volatility weight, base volatility, and volatility transition prior variance).
- two observation parameters: bias and temmperature.
We refer the interested reader to [Mathys et al. 2011](http://www.frontiersin.org/human_neuroscience/10.3389/fnhum.2011.00039/abstract).

The script `demo_dynLearningRate.m` uses the above model to demonstrate how one can use inversion diagnostics to improve existing models, in the aim of progressively capturing the complexity of (hidden) computational mechanisms.