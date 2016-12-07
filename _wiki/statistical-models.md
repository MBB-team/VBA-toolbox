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


# 1D-RFT: Random Field Theory for the multiple comparison problem

Statistical data analyses may require performing multiple tests, e.g., across time samples within a peri-stimulus time window.
In neuroimpaging (e.g., fMRI), the problem of correcting for multiple comparisons across voxels has been solved using random field theory or RFT (Friston et al., 1991, 1994, 1996; Worsley et al., 1992). RFT provides a general expression for the probability of topological features in statistical maps under the null hypothesis, such as the number of peaks above some threshold. It controls for the family-wise error rate or FWER, i.e. the probability of detecting one or more false positives across the entire search volume. This allows one to make valid peak-level and/or cluster-level inferences that account for spatial dependences between voxels, without having to compromise statistical power (as, e.g., a Bonferroni correction would). More precisely, RFT corrects p-values of local peaks in proportion to the estimated roughness of the underlying continuous random field (Kiebel et al., 1999).
VBA includes a simple version of RFT, which obtains when applied to 1D signals (e.g., intra-EEG traces, eyetracking data, skin conductance responses, etc...). It is based upon two main functions:

- **`RFT_main.m`**: this is a generic call to RFT, which can be tailored to any user-specific application

- **`RFT_GLM_contrast.m`**: this applies RFT to GLM-based contrast inference. We will describe an exmaple application below.

Let us assume that our experiment consists in a 2x2 factorial design, with 8 trials per design cell. On each trial, we measure some peri-stimulus response, e.g., a skin conductance response, which has 10^3 time samples. We want to infer on when, in peri-stimulus time, there is a significant interaction of our two experimental factors.
First, the corresponding design matrix and contrast would look something like this:

```matlab
X = kron(eye(4),ones(8,1));
c = [1;-1;1;-1];
```

Let us simulate data under the null (we smooth the noise to take adavntage of RFT's power):

```matlab
L = 1e3; % size of the 1D field
kernel = exp(-0.5.*([1:L]-L/2).^2/(8*2.355)^2); % smoothing kernel (here: FWHM = 8 time samples)
kernel = kernel./sum(kernel);
e = zeros(L,8*4);
for i=1:8*4
   e(:,i) = conv(randn(L,1),kernel,'same'); % smooth residuals
end
b = zeros(4,L); % effect sizes for each time sample (here, no effect)
y = X*b + e';
```

Now let's apply RFT to solve the multiple comparison problem (across time samples):

```matlab
[stat,out] = RFT_GLM_contrast(X,y,c,'t',1,1);
```

In brief, `RFT_GLM_contrast` (i) computes a 1D statistical field composed of Student's t summary statistic for the contrast `c` sampled at each peri-stimulus time point, and (ii) applies RFt to correct for the multiple comparison problem across time samples.
Let us eyeball the ensuing graphical output :

![]({{ site.baseurl }}/images/wiki/1D-RFT.jpg)


> **Upper panel**: The statistical t-field (y-axis) is plotted against time (x-axis). Local peaks are highlighted in red (if the corrected p-value does not reach significance, here: FWER=5%) or in green (if the corrected p-value reaches significance). The same colour-coding applies to upcrossing clusters (for cluster-level inference). **Middle panel**: RFT analysis summary (essentially: expectations, under the null, of features of the sample field). **Lower panel**: list of corrected p-values (set-, cluster- and peak- level inferences). NB: the column "location" relates to local peaks.

All summary statistics are stored in the `out` structure.

> **Tip**: right-liclicking on either local paeks or upcrossing clusters provides a summary of corrected and uncorrected p-values!


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

