---
title: "Multisession"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

In some experiments, a same subject will perform multiple blocs of the same task. It this case, it will often be better to fit all data at once to maximize the precision of the parameter estimates. However, some parameters can differ accross the sessions, like the initial hidden states that should be reinitialized at the beginning of each bloc.

The toolbox offers a simple way to deal with this situation, as demonstrated in the ```demo_multisession``` script.

# Spliting data into sessions

Let's say you have 3 blocs of 40 trials each. You first have to concatenate all your data in order to get one observation vector or matrix with 3x40=120 columns: 

```
y = [y_bloc1 y_bloc2 y_bloc3] ;
```

Then, you simply have to indicate in the ```options``` structure how the toolbox should split the observations by giving how many timepoints are in each session:

```
options.multisession.split = [40 40 40] ;
```

# Fixing parameters across sessions

By default, spliting your data into multiple sessions will duplicate all parameters (evolution $$\theta$$, observation $$\phi$$ and initial hidden state $$X0$$) such that each session has its own set (with identical priors). You can however fix some parameters to keep them constant across all sessions. 

For example, if you have one observation parameter corresponding ie. to the choice temperature you do not expect to vary across blocs:

```
% fix the first observation parameter
 options.multisession.fixed.phi = 1 ; 
```

If you now want to fix the second and third evolution parameters:

```
% fix the first second and third evolution parameter
 options.multisession.fixed.theta = [2 3] ; 
```

Note that by fixing or not pecific parameters you can easily generate different models that can then be compared to test for a session effect.

# Reading the inversion results

In order to estimate your multisession model, just call ```VBA_NLStateSpaceModel``` as usual. The posterior structure will contain all parameter estimates including the fixed parameters and those duplicated for each sessions. To help you reading out these results, the toolbox also output one posterior structure for each session:

```
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

If you fix some parameters, then the estimated value will be the same in all the ```perSession``` structures.