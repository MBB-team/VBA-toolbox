---
title: "Intertemporal choice"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Would you prefer a can of beer today or a bottle of champagne in one week? Intertemporal choices, involving trade-offs between short-term and long-term outcomes, are pervasive in everyday life. The propensity to favor short-term pleasures defines a form of impulsivity that may have dramatic consequences on professional careers or family relationships. How can some people resist the attraction of short-term pleasures and pursue long-term goals, while others easily succumb and compromise their ultimate expectations?
Intertemporal choice models attempt to capture the idiosyncratic temporal discounting of value. Typically, they are static "utility" models, which come in many flavours. Below, we briefly describe the mathematical form of two different variants of delay discounting models:

- Hyperbolic discount:
  \\[U(r,t)=\frac{r}{1+\lambda t} \text{ with } \lambda>0\\]

- Exponential discount:
  \\[U(r,t)=r e^{-\lambda t} \text{ with }\lambda>0\\]

Here, $$U(r,t)$$ is the utility of an item or action that is associated with the reward $$r$$, which will be delivered at time $$t$$.

Such models typically measure peoples' temporal preferences through the discount rate $$\lambda>0$$. Intertemporal choices are modelled as a softmax mapping of the difference in utiltiy of the two alternative options, i.e.: $$U(r,t)-U(r',t')$$. Non-systematic deviations from the model predictions are captured through the behavioural temperature of the softmax mapping.

The script `demo_discounting` exemplifies an analysis of intertemporal choice data that aim at discriminating whether the delay discounting is hyperbolic or exponential.
