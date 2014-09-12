---
title: "Behavioural/Cognitive models"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

The VBA toolbox contains a few behavioural/cognitive models, which can be used to interpret experimental data in the context of perception, learning and/or decision making studies.
Below, we give a few examples of such models. We briefly expose the main theoretical and experimental issues, and point to the relevant functions of the VBA toolbox (demonstration scripts, evolution/observation functions, etc...).

# Reinforcement learning models

These models attempt to capture the main characteristics of operant learning, i.e. how people adapt their behavioural response from environmental feedback (e.g. rewards and punishments).
The central concept in most reinforcement learning models is value, which quantifies how rewarding is a given action or item.
Empirically speaking, one can access these values by many means: subjective verbal report, vegetative responses (e.g., skin conductance or pupil dilation), or decision making. Note that two fundamental aspects of behaviour are driven by value: energy expenditure (one spends more effort when more value is at stake) and explicit choices (one chooses the alternative that has the highest value). Another aspect of value is that it can be learned, e.g. through reinforcement. Reinforcement learning models essentially capture the (possibly changing) action-outcome contingencies that drive behaviour.
A classical example of such models is the so-called "Q-learning" model, whose mathematical details are described [here]({{ site.baseurl }}/wiki/Fast-demo-Q-learning-model)). In its simplest form, its constituents are:

- a set of (action/item) value states. In two-armed bandit problems, there are two of these. In general, there will be as many values as there are available actions. Some behavioural biases can be captured by the initial conditions on these states.
- a learning rate. This (evolution) parameter controls the impact of prediction error onto the value update. Note that one may want o ask whether the learning rate depends upon experimental factors (pathological condition, gain/loss domains, etc...)
- a behavioural temperature. This (observation) parameter controls the amount of noise there is in peoples' choices. Note that such noise is not only related to exploration: model residuals are also captured by the behavioural temperature.

A demonstration script for Q-learning (`demo_Qlearning.m`) is described [here]({{ site.baseurl }}/wiki/Fast-demo-Q-learning-model).

# Bayesian learning models

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

# Intertemporal choice models

Would you prefer a can of beer today or a bottle of champagne in one week? Intertemporal choices, involving trade-offs between short-term and long-term outcomes, are pervasive in everyday life. The propensity to favor short-term pleasures defines a form of impulsivity that may have dramatic consequences on professional careers or family relationships. How can some people resist the attraction of short-term pleasures and pursue long-term goals, while others easily succumb and compromise their ultimate expectations?
Intertemporal choice models attempt to capture the idiosyncratic temporal discounting of value. Typically, they are static "utility" models, which come in many flavours. Below, we briefly describe the mathematical form of two different variants of delay discounting models:

- Hyperbolic discount:
  \\[U(r,t)=\frac{r}{1+\lambda t} \text{ with } \lambda>0\\]

- Exponential discount:
  \\[U(r,t)=r e^{-\lambda t} \text{ with }\lambda>0\\]

Here, $$U(r,t)$$ is the utility of an item or action that is associated with the reward $$r$$, which will be delivered at time $$t$$.

Such models typically measure peoples' temporal preferences through the discount rate $$\lambda>0$$. Intertemporal choices are modelled as a softmax mapping of the difference in utiltiy of the two alternative options, i.e.: $$U(r,t)-U(r',t')$$. Non-systematic deviations from the model predictions are captured through the behavioural temperature of the softmax mapping.

The script `demo_discounting` exemplifies an analysis of intertemporal choice data that aim at discriminating whether the delay discounting is hyperbolic or exponential.
