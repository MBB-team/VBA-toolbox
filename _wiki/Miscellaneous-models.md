---
title: "Miscellaneous models"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

The toolbox also contains additional models that may not be specific to neurobiological and/or behavioural data. These are typically motivated from generic statistical and/or dynamical considerations, which are briefly described below.

# Statistical models

## Classical GLM data analysis

The general linear model (GLM) is a statistical linear model. The GLM incorporates a number of different statistical models: ANOVA, ANCOVA, MANOVA, MANCOVA, ordinary linear regression, t-test and F-test. It is a generalization of multiple linear regression model to the case of more than one dependent variable. The following script demonstrates VBA's GLM classical test functionality:

```matlab
X  = randn(32, 4)  ;
b  = [1; 0; 0; 0 ] ;
y  = X*b ;
y  = y + randn(size(y)) ;
c  = [1; 0; 0; 0 ] ;
pv = GLM_contrast(X, y, c, type) ;
```
where `X` is the design matrix, `y` is the data matrix, `c` is the contrast matrix, `type` is flag for the test (t-test or F-test) and `pv` is the p-value of the corresponding test.

Let us eyeball the graphical output of the function `GLM_contrast.m`:

![]({{ site.baseurl }}/images/wiki/tabs/glm1.jpg)

> **Upper-left panel**: observed data (y-axis) plotted against predicted data (x-axis). NB: The percentage of explained variance is indicated (the adjusted R^2 is the percentage of variance explained by the contrast of interest). **Middle-left panel**: observed and predicted data (y-axis) plotted against data dimensions (x-axis). **Lower-left panel**: parameter's correlation matrix. **Upper-right panel**: parameter estimates plus or minus one standard deviation. NB: Right-clicking on barplots allows one to access the results individual significance F-tests. **Middle-right panel**: contrast of interest. **Lower-right panel**: design matrix. Note that all descriptive statistics can be retrieved from optional output arguments.

## Binary data classification

Strictly speaking, VBa's default generative model copes with continuous data, for which there is a natural distance metric. Now if the data is categorical, there is no such natural metric, and one has to resort to probability distributions dealing with discrete events. For example, binary data can be treated as binomial (Bernouilli) samples, whose sufficient statistic (first-order moment) is given by the observation function. An interesting application is binary data classification. In the linear case, the observation mapping reduces to a linear mixture of explanatory variables (i.e. features) passed through a sigmoid mapping. In this case, VBA's objective is to estimate feature weights that generalize over trials. The script `demo_bin.m` exemplifies such procedure, and demonstrates the relation between classical cross-validation approaches and Bayesian model comparison. See also the demonstration script `demo_logistic.m`.

## Between-conditions group BMS

Until recently, condition and group effects have been addressed by performing random effects BMS independently for the different conditions, and then checking anecdotally to see whether the results of random effects BMS were consistent. This approach is limited, because it does not test the hypothesis that the same model describes the two conditions. VBA provides the inference machinery to evaluate the evidence for a difference – in terms of models – between conditions. This is exemplified in the demosntration script `demo_groupbtw.m`. Note that the main call function for such analysis writes as follows:

```matlab
[ep, out] = VBA-groupBMCbtw(L)
```
where `L` is the log-evidence ND-array, `ep` is the appropriate exceedance probability, and `out` gathers diagnostic metrics.

Readers interested in standard group-BMS are encouraged to run the script `demo_bmc4glm.m`, which relies demonstrates random-effect group BMS in the context of nested models. NB: the model evidence is evaluated at the frequentist limit, using a dedicated function, as follows:

```matlab
[lev] = lev_GLM(y, X)
```
where `y` is the data vector, `X` is the design matrix, and `lev` is the log-model evidence evaluated at the frequentist limit (flat priors).

## Kalman filter/smoother

The script `demo_KalmanSmoother.m` demonstrates the smoothing properties of the Kalman lagged filter.
Here, hidden states follow a triangular wave, whose observation is perturbed with white noise. Critically, we render the inversion scheme blind during half a period of the oscillation. We then invert the model (under AR priors on hidden states), with and without large backward lag.

# Dynamical models

The VBA toolbox includes models that serve as benchmarks for model inversion. Some of these are listed below:

## Van der Pol oscillator

This model captures the essence of action potential dynamics, which essentially consist in stereotypical responses to perturbations, i.e. a stable limit cycle. See `demo_VanDerPol.m`.

![]({{ site.baseurl }}/images/wiki/tabs/VdP1.jpg)

## Lorenz attractor

The Lorenz attractor was originally proposed as a simplified version of the Navier-Stokes equations, in the context of meteorological fluid dynamics. The Lorenz attractor models the autonomous formation of convection cells, whose dynamics are parameterized using three parameters: the Rayleigh number (which characterizes the fluid viscosity),  the Prandtl number (which measures the efficacy of heat transport through the boundary layer) and a dissipative coefficient. When the Rayleigh number is bigger than one, the system has two symmetrical fixed points , which act as a pair of local attractors. For certain parameter values, the Lorenz attractor exhibits chaotic behaviour on a butterfly shaped "strange attractor". For almost any initial conditions (other than the fixed points), the trajectory unfolds on the attractor. The path begins spiralling onto one wing and then jumps to the other and back in a chaotic way. See `demo_Lorenz.m`.

![]({{ site.baseurl }}/images/wiki/tabs/lorenz1.jpg)

## Double-well bistable system

The double-well potential system models a dissipative system, whose potential energy is a quadratic (double-well) function of position. As a consequence, the system is bistable with two basins of attraction to two stable fixed points. In its deterministic variant, the system ends up spiralling around one or the other attractors, depending on its initial conditions and the magnitude of a damping force or dissipative term. Because we consider state-noise, the stochastic DCM can switch (tunnel) from one basin to the other, which leads to itinerant behaviour; this is why the double-well system can be used to model bistable perception. See `demo_doubleWell.m`.

![]({{ site.baseurl }}/images/wiki/tabs/dbw1.jpg)

## Rossler attractor

This is a model of a continuous-time dynamical system that exhibits chaotic dynamics associated with the fractal properties of the attractor. See `demo_Rossler.m`.

![]({{ site.baseurl }}/images/wiki/tabs/rossler1.jpg)

## Henon map

This is one of the simplest discrete-time dynamical systems that exhibit chaotic behavior. See `demo_Henon.m`.

![]({{ site.baseurl }}/images/wiki/tabs/henon1.jpg)
