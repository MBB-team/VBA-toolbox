---
title: "Model inversion in 4 steps"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

This page summarizes the steps required for performing a model inversion with the core routines of the VBA toolbox. Recall that the main model inversion routine is `VBA_NLStateSpaceModel.m`, whose inputs/outputs are defined below. It implements a variational Bayesian approach to the inversion of a very generic class of generative models, namely: "nonlinear state-space models". This class of models is described [here]({{ site.baseurl }}/wiki/Structure-of-VBA's-generative-model).

In brief, in the aim of performing a model-based data analysis using VBA:

- one **needs** to define evolution and observation functions.
- one **can** provide further information about the model and/or its inversion (e.g. priors).

> **TIP:** Many demonstration scripts are provided with the toolbox (e.g., see this [fast demo]({{ site.baseurl }}/wiki/Fast-demo-Q-learning-model)). We suggest you go through some of these to get started.

# Step 1: Defining observation/evolution functions

Generative models are defined in terms of **evolution and observation functions**. The VBA toolbox already contains a library of such functions (for models of behavioural and/or physiological responses). In case one may want to use VBA on one's in-house model, one has to write these evolution/observation functions, in compliance with the following i/o:

```matlab
z = function_name(x_t, P, u_t, in) ;
```

- `x_t` : the vector of hidden states at time `t`
- `P` : the vector of parameters (evolution parameters for the evolution function, observation parameters for the observation function)
- `u_t` : the input (experimenter control variable) at time `t`.
- `in` : may contain any extra relevant information (arbitrary)
- `z`: the predicted state (evolution function) or data (observation function).

> The definition of hidden states ($$x$$), parameters ($$\theta$$ and $$\phi$$), hyperparameters ($$\sigma$$ and $$\alpha$$) and inputs ($$u$$), as well as their role in the model, are given [here]({{ site.baseurl }}/wiki/Structure-of-VBA's-generative-model).

Note that, except the `in` input to evolution/observation functions (which can be anything), all other i/o variables should be **column vectors**. And the dimension of these vectors should correspond to those given in the `dim` structure (see below)!

> **TIP:** In fact, you can also specify the gradients of evolution/observation functions w.r.t states and paraemters as two additional (optional) output variables. This is is useful for **accelerating VBA's inversions** (you can earn up to 2 or 3 orders of magnitude in terms of computational time). For example:
```matlab
function [g,dgdx,dgdP] = g_dummy(x,P,u,in)
% dummy observation function with one state and one parameter 
g = x.*P.*u; % "minimal" output
dgdx = P.*u; % gradient wrt states (x)
dgdP = x.*u; % gradient wrt parameters (P)
end
```
Thus, the second and third outputs of evolution/observation functions are optional but, if specified, they should only be used for gradients!


# Step 2 : Setting model inversion options

The VBA model inversion requires the user to specify some additional information:

- **model dimensions** : VBA stores these in the structure `dim`, which contains the following fields:
  - `n` : number of hidden states
  - `p` : output (data) dimension, ie. number of obervations per time sample
  - `n_theta` : number of evolution parameters
  - `n_phi` : number of observation parameters
  - `n_t` : number of time samples

For example, setting:

```matlab
dim.n       = 1 ; % number of hidden states
dim.n_theta = 2 ; % number of evolution parameters
dim.n_phi   = 3 ; % number of observation parameters
```
tells VBA that there are 1 hidden state, 2 evolution parameters and 3 observation parameters.

> **TIP:** The `dim` structure has to match the dimensions of the inputs (more precisely: hidden states and parameters) to the evolution and observation functions. In case you don't use VBA's default priors, they also have to match the dimensions of your priors...

- Other **options**

VBA allows the user to control the inversion using an `options` input structure, which is passed to `VBA_NLStateSpaceModel.m`. These options include, but are not limited to: informing VBA about categorical and/or missing data, setting "micro-time resolution", passing optional arguments to evolution and/or observation functions, etc... [This page]({{ site.baseurl }}/wiki/Controlling-the-inversion-using-VBA-options) provides an (almost) exhaustive list of these options.


# Step 3 : Defining priors

In addition to the evolution and observation functions, specifying the generative model requires the definition of **prior probability distributions** over model unknown variables. These are summarized by sufficient statistics (e.g., mean and variance), which are stored as a matlab structure `priors` that is itself apended to the `options` structure:

- **Observation parameters**
  - `priors.muPhi`: prior mean on $$\phi$$
  - `priors.SigmaPhi`: prior covariance on $$\phi$$
- **Evolution parameters** (only for dynamical systems)
  - `priors.muTheta`: prior mean on $$\theta$$
  - `priors.SigmaTheta`: prior covariance on $$\theta$$
- **Initial conditions** (only for dynamical systems)
  - `priors.muX0`: prior mean on $$x_0$$
  - `priors.SigmaX0`: prior covariance on $$x_0$$
- **Measurement noise precision** (only for continuous data)
  - `priors.a_sigma`: prior shape for the measurement noise precision $$\sigma$$
  - `priors.b_sigma`: prior rate for the measurement noise precision $$\sigma$$
- **State noise precision** (only for dynamical systems)
  - `priors.a_alpha`: prior shape for the state noise precision $$\alpha$$
  - `priors.b_alpha`: prior rate for the state noise precision $$\alpha$$
- **Noise precision matrices**
  - `priors.iQy`: prior precision matrices for the measurement noise (only for gaussian data sources)
  - `priors.iQx`: prior precision matrices for the state noise (only for stochastic systems)

If left unspecified, the `priors` structure is filled in with defaults (typically, i.i.d. zero-mean and unit-variance Gaussian densities, except for $$\sigma$$ and $$\alpha$$). For example, setting:

```matlab
priors.muPhi = zeros(dim.n_phi,1); % prior mean (obs params)
priors.SigmaPhi = eye(dim.n_phi); % prior covariance (obs params)
```
effectively defines a $$N\left( 0,I \right)$$ i.i.d. (zero-mean, unit-variance) normal density on observation parameters.

> **Tip:** By default, VBA assumes that dynamical systems are deterministic. This is done by forcing an infinite prior precision on state noise. However, VBA identification of  **stochastic** systems can be performed by setting a finite state noise precision.
>
> For example, setting:
>
>```matlab
priors.a_alpha = 1; % prior shape param (state noise precision) 
priors.b_alpha = 1; % prior rate param (state noise precision)
```
effectively assumes that state noise precision $$\alpha$$ is a unit-mean and unit-variance [Gamma variable](https://en.wikipedia.org/wiki/Gamma_distribution). This allows non-zero state noise to enter and perturb the hidden states' dynamics...

One then fills in the `priors` field of the `options` structure, as follows:

```matlab
options.priors = priors ; % store 'priors' in VBA's 'options' input structure
```

> In any bayesian data analysis, setting the priors is a subtle issue. From a classical perspective, priors induce a systematic bias in parameter estimation (but remember the "[bias-variance trade-off](https://en.wikipedia.org/wiki/Bias%E2%80%93variance_tradeoff)" of machine learning). More importantly, when data is of insufficient quantity and/or quality, priors partly influence bayesian model comparison, which is based upon the [marginal likelihood](https://en.wikipedia.org/wiki/Marginal_likelihood). At this point, suffices to say that VBA enables so-called ["empirical Bayes"](https://en.wikipedia.org/wiki/Empirical_Bayes_method) approaches, in which the priors are estimated from the data. This is explained on [this page]({{ site.baseurl }}/wiki/VBA-MFX)...


# Step 4 : Inverting the model

Having completed steps 1 to 3, one simply calls the main **VB model inversion** routine, namely `VBA_NLStateSpaceModel.m`, as follows:

```matlab
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options)
```

Its input arguments are:

- the data `y`

- the input `u` (can be left empty)

- `f_fname`: the name/handle of the evolution function (left empty for static models)

- `g_fname`: the name/handle of the observation function

- `dim`: the dimensions of the model variables

- `options` (can be left empty)

Its output arguments are:

- `posterior`: a matlab structure that has the same format as the `priors` above (i.e. stores the first- and second-order moments of the posterior densities of all unknown variables in the model)

- `out`: a matlab structure that summarizes some diagnostics of the model VBA inversion. NB: the (lower bound to) the model evidence is stored in `out.F`.

> A simple and fast demonstration of a VBA model inversion (for a Q-learning model) is described [here]({{ site.baseurl }}/wiki/Fast-demo-Q-learning-model). A description of VBA's graphical output is described [here]({{ site.baseurl }}/wiki/VBA-graphical-output).

