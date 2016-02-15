---
title: "Statistical models"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

# Classical GLM data analysis

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

**Tip**: defining contrasts for, e.g., main effects of multi-level factors can be tedious. However, VBA has an in-built function for doing exactly this: `Contrast_MEbins.m`. It outputs the contrast matrix corresponding to an F-test of the main effect of a given experimental factor with n levels. 

# GLM with missing data

Missing data occur when the value of either the dependent or the independent variables (or both) are unknown or lost. This is an issue because the GLM cannot be fitted using the usual procedure (cf. missing values in the design matrix). VBA relies on the variational Bayesian algorithm to estimate likely values for missing data in addition to fitting the GLM parameters. See the demo: `demo_GLM_missingData.m`.


# Mediation analysis

One may be willing to address questions about how key variables mediate the impact of (controlled) stimuli onto (measured) experimental outcomes. For example, one may want to assess the impact of brain activity onto behavioural outcomes above and beyond the effect of psychophysical manipulations. Statistical models of mediated effects place the mediator variable (M) at the interplay between the independent variable (X) and the dependent variable (Y).
On the one hand, VBA proposes to approach such mediation analyses from a Bayesian perspective. The latter perspective reduces to a bayesian model comparison, whereby one evaluates the evidence in favour of a model that includes a serial chain of causality (X -> M -> Y) against a "common cause" model (X -> M and X-> Y).
On the other hand, VBA comprises self-contained routines for performing classical tests of mediation (e.g., Sobel tests). In brief, this consists in testsing for the significance of the product of path coefficients in the serial model above. Such test can be performed using the function `mediationAnalysis0.m` (see the demo: `demo_mediation.m`.


# Binary data classification

Strictly speaking, VBa's default generative model copes with continuous data, for which there is a natural distance metric. Now if the data is categorical, there is no such natural metric, and one has to resort to probability distributions dealing with discrete events. For example, binary data can be treated as binomial (Bernouilli) samples, whose sufficient statistic (first-order moment) is given by the observation function. An interesting application is binary data classification. In the linear case, the observation mapping reduces to a linear mixture of explanatory variables (i.e. features) passed through a sigmoid mapping. In this case, VBA's objective is to estimate feature weights that generalize over trials. The script `demo_bin.m` exemplifies such procedure, and demonstrates the relation between classical cross-validation approaches and Bayesian model comparison. See also the demonstration script `demo_logistic.m`.

# Kalman filter/smoother

The script `demo_KalmanSmoother.m` demonstrates the smoothing properties of the Kalman lagged filter.
Here, hidden states follow a triangular wave, whose observation is perturbed with white noise. Critically, we render the inversion scheme blind during half a period of the oscillation. We then invert the model (under AR priors on hidden states), with and without large backward lag.

# Probabilistic clustering

Any Bayesian data analysis relies upon a generative model, i.e. a probabilistic description of the mechanisms by which observed data are generated. The main VBA functions deal with a certain class of generative models, namey: state-space models with unknown evolution/observation parameters. But it also offers a set of VB routines developped to invert mixtures of Gaussian -MoG- and of binomial -MoB- densities, which are described below. These can be used to perform data clustering (i.e., blind data classification). 
The '\classification' folder contains the functions that deal with MoG and MoB data (see below).

The MoG model has been widely used in the machine learning and statistics community. It assumes that the observed data y is actually partitioned into K classes or subsets of points, each one of which can be well described with a multivariate Gaussian density. More formally, the MoG generative model assumes that n binomial vectors are sampled from a Dirichlet density. These binomial vectors have dimension K, and are such that their non-zero entry indicates which class each point belongs to. These so-called labels then switch on and off a set of K multivariate Gaussian densities with different first- and second- order moments, from which the observed data themselves are sampled. If these densities are sufficiently different, then the overall datasets enjoys a clustering structure, when projected onto an appropriate subspace.
NB: Any arbitrary density over continuous data can be described in terms of a MoG, given a sufficient number of components or classes.
The function `VBA_MoG.m` inverts, using a variational Bayesian scheme, the MoG model. Given the data y, and the maximum number of classes K, it estimates the most plausible number of classes, the labels, and the first- and second- order moments of each Gaussian density. In addition, it also returns the model evidence, which can be useful for model comparison purposes (see the demo: `demo_GMM.m`).

Now suppose you design an experiment, in which n subjects are asked a number yes/no questions from a questionnaire. Suppose these subjects are grouped into K categories, which are defined in terms of how likely its member are to answer 'yes' to each of the questions. You have just defined a MoB model, which you can use to disclose the categories of subjects from their profile of responses to the questionnaire. This is what the function '
`MixtureOfBinomials.m` does, using either a Gibbs sampling algorithm or a variational Bayesian algorithm. In both cases, the function returns the data labels, the first- order moment profile for each category and the model evidence. NB: in the MoB case, the inversion schemes cannot automatically eliminate the unnecessary classes in the model. Therefore, (Bayesian) model comparison is mandatory to estimate the number of classes in the data, if it is unknown to the user. See the demo: `demo_BMM.m`.

