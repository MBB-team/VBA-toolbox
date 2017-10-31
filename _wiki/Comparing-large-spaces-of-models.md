---
title: "Bayesian model selection with large model spaces"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Bayesian model selection (BMS) allows one to compare any set of models, provided the experimental data is the same. In particular, one may want to compare models that differ in terms of whether or not a given parameter is zero. One would say that these models are "[nested](https://en.wikipedia.org/wiki/Statistical_model#Nested_models)", i.e. the "reduced" model obtains from the "full" model by setting the prior mean and variance of the corresponding parameter to zero. This is the BMS equivalent to classical significance testing (c.f. F-test), where the "null" is simply the "reduced" model.

Interestingly, one can show that the Bayes factor of nested models can be simply derived from the prior and posterior densities of the "full" model, as follows:

$$\frac{p(y\mid H_0)}{p(y\mid H_1)}=\frac{p(\vartheta=0\mid y,H_1)}{p(\vartheta=0\mid H_1)}$$

where $$H_1$$ is the "full" model, and $$H_0$$ is the "reduced" model (i.e. with the parameter $$\vartheta$$ fixed to zero). The right-hand side of the equation above is known as the **Savage-Dickey ratio**.

This implies that one only needs to invert the full model to test the significance of its parameters. The following script demonstrates the approach:

```matlab
% invert full model
[posterior, out] = VBA_NLStateSpaceModel(y, u, f_fname, g_fname, dim, options) ;
   
% define reduced model
rprior = out.options.priors ;
rprior.muPhi(1)    = 0 ;
rprior.SigmaPhi(1) = 0 ;

% Compute reduced model evidence
[Fr, rpost] = VBA_SavageDickey(posterior, out.options.priors, out.F, dim, rprior) ;
```

The first line is VBA's inversion of the full model (arbitrary data and model structure). The next lines create the reduced model's priors from the full model's, by zeroing the prior mean and variance of the first observation parameter. Then, both the log-evidence of the reduced model (`Fr`) and its posterior density (`rpost`) are derived. Note that posterior correlation among parameters under the full model induce changes in the marginal posterior densities of the non-zero parameters of the reduced model.

The key insight here is that the use of Savage-Dickey ratios is orders of magnitude faster than a proper model inversion. This means that one can use this to compare large spaces of nested models, by looping over combinations of model parameters. However, Savage-Dickey ratios only provide approximated inference!

> For further details, please refer to, e.g., [Rosa et al. (2013)](https://www.ncbi.nlm.nih.gov/pubmed/21459150)

