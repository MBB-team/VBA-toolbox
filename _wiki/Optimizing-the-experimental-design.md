---
title: "Experimental design optimization"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

This section describes how to optimize the experimental design with the aim of either estimating model parameters or comparing models. First, we recall the definition of design optimality scores. Then, we suggest on-line extensions to the approach for adapative design strategies.

# Design optimality scores

Optimizing the design in the context of, e.g., experimental psychology studies, amounts to identifying the subset of conditions and stimuli ($$u$$) that yields the highest statistical power, under a variety of practical constraints. This requires being able to predict experimental data under models' assumptions, including potential confounds that may mask the effects of interest. Design optimization can become a difficult issue whenever the impact of experimental factors onto measurements (through the model) is non-trivial and/or uncertain (cf. unknown model parameters). This motivates the use of automatic design optimization.

The VBA toolbox can handle two classes of problems, namely optimizing the system's input $$u$$ with respect to either parameter estimation or model selection. These two problems correspond to two different objectives, which can be formalized in terms of statistical loss functions:

## Parameter estimation

How should one set the inputs $$u$$, such that measured experimental data eventually yield accurate parameter estimates? One first has to define "estimation accuracy". Let $$\hat{\theta} = E \big[ \theta\mid y \big]$$ be the posterior estimate of unknown model parmaeters $$\theta$$, given experimental data $$y$$. The expected estimation error is given by:

$$ E \big[( \hat{\theta} -\theta)^2 \mid y \big] = V \big[ \theta\mid y \big]$$

where $$V \big[ \theta\mid y \big]$$ is the posterir variance over the unknown model parameter $$\theta$$. Thus, maximzing estimation accuracy reduces to minimizing a posteriori variance. This intuition generalizes to models with more than one parameter. In this case, one usually minimizes the trace of the expected posterior matrix (cf. so-called "A-optimality"). Note that, for non-linear generative models, the posterior variance will depend upon the (not yet observed) data $$y$$. One then needs to derive the expected posterior variance, where the expectation is taken under the marginal distribution of the data. This can be done by evaluating the design efficiency for each design, as follows:

```matlab
e = VBA_designEfficiency(f_fname,g_fname,dim,options,u,'parameters')
```
where `u` is the time series of experimental control variables (inputs) that defines the design and `e` is the design efficiency (minus the trace of the expected posterior covariance matrix).

## Model selection

In this context, we argue that one should choose among experimental designs according to their induced model selection error rate. This can be done by choosing the input $$u$$ that minimizes the so-called "Laplace-Chernoff risk"  (cf. [Daunizeau et al. 2011](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002280)) $$b_LC$$ , which yields a lower bound on the model selection error rate:

$$b_{LC}(u) = 1-\frac{1}{2}\log\:\left(\frac{\Delta g(u)^2}{4\tilde Q(u)}+1\right) \quad \text{if}\; \tilde Q_1(u) \approx \tilde Q_2(u) \equiv \tilde Q_(u) $$

![]({{ site.baseurl }}/images/wiki/optim/optim0.jpg)

Here, prior predictive densities (y-axis) of two different models are plotted over possible data values (x-axis; unidimensional data). The expected model selection error rate $$p(e=1|u)$$ increases with the similarity of the two probabilistic model predictions. In fact, the Laplace-Chernoff risk $$b_LC$$ measures the statistical similarity of the prior predictive densities, as a function of their 1st- and 2nd-order moments.

The numerical derivation of the design efficiency for model comparison can be done as follows:

```matlab
[e,out] = VBA_designEfficiency(f_fname,g_fname,dim,options,u,'models')
```
where `f_fname`, `g_fname`, `dim`, `options` are nx1 cell arrays (n models), `e` is design efficiency (i.e. minus the Chernoff risk) and `out` is a structure containing diagnostic variables (e.g.: upper bound on selection error probability, 1st- and 2nd-order moments of the Laplace approximation to the prior predictive density..).

One can then either compare different designs on the basis of their efficiency (e.g. on a predefined set of inputs $$u$$), or perform numerical optimization of the design efficiency w.r.t. $$u$$ (or some parametric form of it).


## On-line adaptive designs

In the context of experimental psychophysics, adaptive designs such as "stair-case" methods are used to, e.g., estimate some individual sensory detection or discrimination threshold. Such procedures operate in real-time in the sense that the next stimulation depends on the previous behavioral response and is computed in order to optimize model fitting. More generally, adaptive (on-line) designs can be used to improve on three problems: (i) model parameter estimation; (ii) hypothesis testing (or BMS); (iii) choosing the duration of the experiment (e.g., the number of trials).

We have implemented an example of on-line adaptive design in the demonstration script `demo_binomial_AdaptDesign.m`. This demo simulates a psychophysics paradigm similar to a signal detection task, whereby the detection probability is a sigmoidal function of the stimulus contrast (which is the design control variable). The gol of the experiment is to estimate the sigmoid inflexion point (detection threshold) and the sigmoid slope (d-prime). Here, the design is adapted online, in the aim of providing the most efficient estimate of these model parameters, given trial-by-trial subjects' binary choice data (`y=1`: "seen", `y=0`: "unseen").

![]({{ site.baseurl }}/images/wiki/optim/optim1.jpg)

On the upper-left graph, one can see the evolution of the posterior credible intervals over the model parameters (y-axis; blue: sigmoid slope, green: inflexion point) as a function of trial index (x-axis). These converge very quickly (and precisely) around the simulated values. The upper-right graph shows the estimated response curve, in terms of the detection probability (y-axis) as a function of stimulus contrast (x-axis). The lower-left graph plots the design efficiency (y-axis) as a function of stimulus contrast (x-axis), at the last trial. One can see that the design is most efficient when the sigmoid function is sampled around the maximal curvature regions. The lower-right graph depicts the histogram of design control variables over trials, which has essentially focused the sigmoid sampling around the inflexion points.

Note that the experiment could have been stopped much before (after about 20 trials), on the basis of the convergence of the design efficiency score:

![]({{ site.baseurl }}/images/wiki/optim/optim2.jpg)


Early-stopping rules can be enforced by setting a threshold on design efficiency, in the aim of performing short and efficient experiments...


