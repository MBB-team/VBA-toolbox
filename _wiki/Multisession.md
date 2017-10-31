---
title: "Multisession inversions"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

In some experiments, a given subject may perform multiple blocks of the same task. It this case, one may want to fit all data at once to maximize the precision of the parameter estimates. However, some parameters can differ accross the sessions, such as the initial hidden states (these should be reinitialized at the beginning of each block). More generally, one may want to let some parameters vary over sessions, while others are kept constant (e.g., to test for a specific session effect).

The toolbox offers a simple way to deal with situations of this sort, as demonstrated in the ```demo_multisession``` script.

# Spliting data into sessions

Let's say you have 3 blocs of 40 trials each. You first have to concatenate all your data in order to get one observation vector or matrix with 3x40=120 columns: 

```matlab
y = [y_block1 y_block2 y_block3] ;
```

Then, you simply have to indicate in VBA's ```options``` structure how the observations should be patitionned, by specifying the number of time samples there is for each block:

```matlab
options.multisession.split = [40 40 40] ;
```

# Fixing parameters across sessions

By default, spliting your data into multiple sessions will duplicate all parameters (evolution $$\theta$$, observation $$\phi$$ and initial hidden state $$x_0$$) such that each session has its own set (with identical priors). In this case, there is no assimilation of information from one session (or block) to the next. However, you can tell VBA to fix some parameters, i.e. to keep them constant across all sessions. In turn, these parameters will be estimated given all the available sessions.

For example, let's say that the first observation parameter (corresponding, e.g., to the choice temperature) is not expected to vary across sessions:

```matlab
% fix the first observation parameter
options.multisession.fixed.phi = 1 ; 
```

If you now want to fix the second and third evolution parameters:

```matlab
% fix the first second and third evolution parameters
options.multisession.fixed.theta = [2 3] ; 
```

Note that by fixing or not specific parameters you can easily generate different models that can then be compared to test for a session effect.

# Reading the inversion results

In order to estimate your multisession model, just call ```VBA_NLStateSpaceModel``` as usual, but with the `options` structure set as above. The `posterior` structure will contain the usual outputs, plus session-specific posterior structures:

```matlab
% one structure for each session 
 posterior.perSession 
 ans = 
 1x3 struct array with fields:
    muX0
    SigmaX0
    muTheta
    SigmaTheta
    muPhi
    SigmaPhi
    a_sigma
    b_sigma
    a_alpha
    b_alpha
    iQy
    iQx
    muX
    SigmaX
```

If you fix some parameters, then the estimated value will be the same in all the ```posterior.perSession``` structures.
