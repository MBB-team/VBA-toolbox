---
title: "Model inversion in 4 steps"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

This page summarizes the steps required for performing a model inversion with the core routines of the VBA toolbox. Recall that the main model inversion routine is `VBA_NLStateSpaceModel.m`, whose inputs/outputs are defined below. It implements a variational Bayesian approach to the inversion of a very generic class of generative models, namely: "nonlinear state-space models". This class of models is described [here]({{ site.baseurl }}/wiki/Structure-of-VBA's-generative-model).

In brief, in the aim of performing a model-based data analysis using VBA:

- one **needs** to define evolution and observation functions, as well as creating the `dim` matlab structure (see below).
- one **can** provide further information about the model and/or its inversion (e.g. priors).

> **TIP:** Many demonstration scripts are provided with the toolbox (e.g., see this [fast demo]({{ site.baseurl }}/wiki/Fast-demo-Q-learning-model)). We suggest you to go through some of these to get started.

# Step 1: Defining observation/evolution functions

Generative models are defined in terms of **evolution and observation functions**. The VBA toolbox already contains a library of such functions (for models of behavioural and/or physiological responses). In case one may want to use VBA on one's in-house model, one has to write these evolution/observation functions, in compliance with the following I/O:

```matlab
z = function_name(x_t, P, u_t, in) ;
```

- `x_t` : the vector of hidden states at time `t`
- `P` : the vector of parameters (evolution parameters for the evolution function, observation parameters for the observation function)
- `u_t` : the input (experimenter control variable) at time `t`.
- `in` : may contain any extra relevant information (arbitrary)
- `z`: the predicted state (evolution function) or data (observation function).

> The definition of hidden states ($$x$$), parameters ($$\theta$$ and $$\phi$$), hyperparameters ($$\sigma$$ and $$\alpha$$) and inputs ($$u$$), as well as their role in the model, are given [here]({{ site.baseurl }}/wiki/Structure-of-VBA's-generative-model).

# Step 2 : Setting model inversion options

The VBA model inversion requires the user to specify some additional information:

- **model dimensions** : `dim`
  - `n` : number of hidden states
  - `p` : output (data) dimension, ie. number of obervations per time sample
  - `n_theta` : number of evolution parameters
  - `n_phi` : number of observation parameters
  - `n_t` : number of time samples

For example, setting:

```matlab
dim.n       = 1 ;
dim.n_theta = 2 ;
dim.n_phi   = 3 ;
```
tells VBA that there are 1 hidden state, 2 evolution parameters and 3 observation parameters.

> **TIP:** other dimensions (`dim.p` and `dim.n_t`) are optional.

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


If left unspecified, the `priors` structure is filled in with defaults (typically, i.i.d. zero-mean and unit-variance Gaussian densities, except for $$\sigma$$ and $$\alpha$$). For example, setting:

```matlab
priors.muPhi    = zeros(dim.n_phi,1) ;
priors.SigmaPhi = eye(dim.n_phi)     ;
```
effectively defines a $$N\left( 0,I \right)$$ i.i.d. (zero-mean, unit-variance) normal density on observation parameters.

> **Tip:** By default, VBA assumes that dynamical systems are deterministic. This is done by forcing an infinite prior precision on state noise. However, identification of  **stochastic** systems can be performed by setting a finite state noise precision.
>
> For example, setting:
>
>```matlab
priors.a_alpha = 1;
priors.b_alpha = 1;
```
effectively assumes that state noise precision $$\alpha$$ is a unit-mean and unit-variance [Gamma variable](https://en.wikipedia.org/wiki/Gamma_distribution).

One then fills in the `priors` field of the `options` structure, as follows:

```matlab
options.priors = priors ;
```


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

