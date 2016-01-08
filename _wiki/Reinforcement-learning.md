---
title: "Reinforcement learning"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

These models attempt to capture the main characteristics of operant learning, i.e. how people adapt their behavioural response from environmental feedback (e.g. rewards and punishments).
The central concept in most reinforcement learning models is value, which quantifies how rewarding is a given action or item.

Empirically speaking, one can access these values by many means: subjective verbal report, vegetative responses (e.g., skin conductance or pupil dilation), or decision making. Note that two fundamental aspects of behaviour are driven by value: energy expenditure (one spends more effort when more value is at stake) and explicit choices (one chooses the alternative that has the highest value). Another aspect of value is that it can be learned, e.g. through reinforcement. Reinforcement learning models essentially capture the (possibly changing) action-outcome contingencies that drive behaviour.

A classical example of such models is the so-called "Q-learning" model, whose mathematical details are described [here]({{ site.baseurl }}/wiki/Fast-demo-Q-learning-model)). In its simplest form, its constituents are:

- a set of (action/item) value states. In two-armed bandit problems, there are two of these. In general, there will be as many values as there are available actions. Some behavioural biases can be captured by the initial conditions on these states.
- a learning rate. This (evolution) parameter controls the impact of prediction error onto the value update. Note that one may want o ask whether the learning rate depends upon experimental factors (pathological condition, gain/loss domains, etc...)
- a behavioural temperature. This (observation) parameter controls the amount of noise there is in peoples' choices. Note that such noise is not only related to exploration: model residuals are also captured by the behavioural temperature.

A demonstration script for Q-learning (`demo_Qlearning.m`) is described [here]({{ site.baseurl }}/wiki/Fast-demo-Q-learning-model).