---
title: "Setting the priors through Empirical Bayes"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


[Empirical Bayes methods](https://en.wikipedia.org/wiki/Empirical_Bayes_method) are procedures for statistical inference in which the prior distribution is estimated from the data. This approach stands in contrast to standard Bayesian methods, for which the prior distribution is fixed before any data are observed. In VBA, the empirical Bayes approach is a fully Bayesian treatment of a [hierarchical model](https://en.wikipedia.org/wiki/Bayesian_hierarchical_modeling) wherein the parameters at the highest level of the hierarchy are [summary statistics](https://en.wikipedia.org/wiki/Summary_statistics) of the group, which are unknown but eventually constrain the likely range of subject-level parameters. 

In other terms, VBA's empirical Bayes approach is formally identical to [mixed-effect](https://en.wikipedia.org/wiki/Mixed_model) modelling (MFX), whereby the priors on parameters are hierarchically defined, such that subject-level parameters are assumed to be [sampled](https://en.wikipedia.org/wiki/Sample_(statistics)) from a Gaussian density whose mean and variance are unknown but match the moments of the group distribution.

The pseudo-code of VBA's ensuing hierarchical inversion is given below:

```matlab
[LOOP UNTIL CONVERGENCE]
1) for i=1:n (loop over subjects)
      define within-subject priors from group-level summary statistics
      perform within-subject model inversion
  end
2) update group-level summary statistics from posterior within-subject summary statistics
```

Over the algorithm iterations, within-subject priors are refined and matched to the inferred parent population distribution. Empirical Bayes procedures of this sort learn from group statistics, and thus inform within-subject inversions with each other results. This eventually [shrinks](https://en.wikipedia.org/wiki/Shrinkage_estimator) the within-subject posterior estimate around the estimated group mean...

Of course, there is no need to write specific functions, and VBA performs this analysis automatically. One simply calls the function `VBA_MFX.m`, as follows:

```
[p_sub,o_sub,p_group,o_group] = VBA_MFX(y,u,f_fname,g_fname,dim,options,priors_group)
```

Its intputs and outputs are described below:

- `y`: nsx1 cell array of observations, where ns is the number of subjects in the group
- `u`:  nsx1 cell array of inputs
- `f_fname`/`g_fname`: evolution/observation function handles
- `dim`: structure containing the model dimensions.
- `options`: nsx1 cell array of options structure.
- `priors_group`: structure containing the prior sufficient statistics on the moments of the parent population distributions (for observation and evolution parameters, as well as for initial conditions, if applicable). See `p_group` subfields below.
- `p_sub`/`o_sub`: nsx1 cell arrays containng the VBA outputs of the within-subject model inversions.
- `p_group`: structure containing the sufficient statistics of the posterior over the moments of the parent population distribution.
- `o_group`: output structure of the VBA_MFX approach. In particular, it contains the Free Energy of the MFX model (for model comparison purposes).

A complete mathematical description of the approach is available [on this page](https://sites.google.com/site/jeandaunizeauswebsite/links/resources).
