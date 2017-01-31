---
title: "Bi-dimensional decisions"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


Subjective [traits](https://en.wikipedia.org/wiki/Trait_theory) or [attitudes](https://en.wikipedia.org/wiki/Attitude_(psychology)) such as "prudence", "impatience" or "laziness" are key determinants of goal-directed behaviour. This is because they determine how people arbitrate between canonical but conflicting decision dimensions, e.g., the prospect of [reward](https://en.wikipedia.org/wiki/Reward_system) and costs such as [risk](https://en.wikipedia.org/wiki/Risk), [delay](https://en.wikipedia.org/wiki/Temporal_discounting) or [effort](https://en.wikipedia.org/wiki/Principle_of_least_effort). For example, high risk devaluation is the hallmark of "prudence", "impatience" is associated with strong delay discounting and "lazy" people find potential rewards not worth the effort.

The common notion behind these types of behaviour is that they result from a cost-benefit trade-off, i.e. a form of [bi-dimensional decision](https://en.wikipedia.org/wiki/Multiple-criteria_decision_analysis). Modelling such decisions is a difficult issue, and VBA is equipped with a few plug-and-play tools...



## The case of inter-temporal choice models

Would you prefer a can of beer today or a bottle of champagne in one week? [Intertemporal choices](https://en.wikipedia.org/wiki/Intertemporal_choice), involving trade-offs between short-term and long-term outcomes, are pervasive in everyday life. The propensity to favor short-term pleasures defines a form of [impulsivity](https://en.wikipedia.org/wiki/Impulsivity) that has been shown to be be a critical determinant of IQ, drug abuse or eating behaviour [(Mischel et al., 1989)](https://www.ncbi.nlm.nih.gov/pubmed/2658056). How can some people resist the attraction of short-term pleasures and pursue long-term goals, while others easily succumb and compromise their ultimate expectations?

Intertemporal choice models attempt to capture the idiosyncratic temporal discounting of value. Typically, they are static "[utility](https://en.wikipedia.org/wiki/Utility)" models, which come in many flavours. Below, we briefly describe the mathematical form of two different variants of delay discounting models:

- Hyperbolic discount:
  \\[U(r,t)=\frac{r}{1+\lambda t} \text{ with } \lambda>0\\]

- Exponential discount:
  \\[U(r,t)=r e^{-\lambda t} \text{ with }\lambda>0\\]

Here, $$U(r,t)$$ is the utility of an item or action that is associated with the reward $$r$$, which will be delivered at time $$t$$.

Such models typically measure peoples' [time preferences](https://en.wikipedia.org/wiki/Time_preference) through the discount rate $$\lambda>0$$. Intertemporal choices are modelled as a [softmax mapping](https://en.wikipedia.org/wiki/Softmax_function) of the difference in utiltiy of the two alternative options, i.e.: $$U(r,t)-U(r',t')$$. Non-systematic deviations from the model predictions are captured through the behavioural temperature of the softmax mapping.

The script `demo_discounting` exemplifies an analysis of intertemporal choice data that aim at discriminating whether the delay discounting is hyperbolic or exponential.


## Agnostic modelling of bi-dimensional decisions

Other forms of cost-benefit arbitrages have been proposed for, e.g., attitudes toward risk and effort. But the typical issue in bi-dimensional decisions is that there is no simple way of visualizing the utility profile, if any. All we have is a series of choices, i.e. some information about whether people prefer this bundle of features over this other one. This is where one can resort to good old physics tricks, in particular: [Fourier decomposition](https://en.wikipedia.org/wiki/Fourier_series). Fourier series provide an approximation of any signal, with a given (limited) resolution.

Let $$U(x,y)$$ be the utility function defined over a 2D domain composed of two orthogonal choice features (e.g., reward and effort). Then, $$U(x,y)$$ can be decomposed onto a Fourier basis function set, as follows:

$$U(x,y) = a_0 + \sum_n{ \sum_m{ (a_{nm} f_{nm}(x,y)+b_{nm} g_{nm}(x,y)+c_{nm} h_{nm}(x,y)+d_{nm} i_{nm}(x,y)) }}$$

where $$(L1,L2)$$ is the size of the domain and the 2D Fourier basis functions are given by:

$$f_{nm}(x,y) = \cos(2n \pi x/L1) \cos(2m \pi y/L2) $$
$$g_{nm}(x,y) = \cos(2n \pi x/L1) \sin(2m \pi y/L2) $$
$$h_{nm}(x,y) = \sin(2n \pi x/L1) \cos(2m \pi y/L2) $$
$$i_{nm}(x,y) = \sin(2n \pi x/L1) \sin(2m \pi y/L2) $$

Here, $$a_{nm}$$, $$b_{nm}$$, $$c_{nm}$$ and $$d_{nm}$$ are unknown Fourier coefficients, to be adjusted to the signal. Here again, this can be done by modelling each choice as a softmax mapping of the difference in utilty of the two alternative options, i.e.: $$U(x,y)-U(x',y')$$, where the utility function has now been replaced by its Fourier decomposition. Then, one can reconstruct the estimated utility profile by inserting the estimated Fourier coefficients in the series above, and display it.

VBA can be used to derive the Fourier basis function set (cf. function `Fourier2DBF.m`) and then  provide estimates of the Fourier coefficients. The script `demo_2DChoices` exemplifies this approach.

Notes:

- the order $$(n,m)$$ of the series controls the effective [resolution](https://en.wikipedia.org/wiki/Image_resolution) of the reconstruction.
- The so-called "[discrete-cosine-transform](https://en.wikipedia.org/wiki/Discrete_cosine_transform)" or DCT may be preferred to the above Fourier basis function set (the function `Fourier2DBF.m` also deals with DCT).










