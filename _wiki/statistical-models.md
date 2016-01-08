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

# Binary data classification

Strictly speaking, VBa's default generative model copes with continuous data, for which there is a natural distance metric. Now if the data is categorical, there is no such natural metric, and one has to resort to probability distributions dealing with discrete events. For example, binary data can be treated as binomial (Bernouilli) samples, whose sufficient statistic (first-order moment) is given by the observation function. An interesting application is binary data classification. In the linear case, the observation mapping reduces to a linear mixture of explanatory variables (i.e. features) passed through a sigmoid mapping. In this case, VBA's objective is to estimate feature weights that generalize over trials. The script `demo_bin.m` exemplifies such procedure, and demonstrates the relation between classical cross-validation approaches and Bayesian model comparison. See also the demonstration script `demo_logistic.m`.

# Kalman filter/smoother

The script `demo_KalmanSmoother.m` demonstrates the smoothing properties of the Kalman lagged filter.
Here, hidden states follow a triangular wave, whose observation is perturbed with white noise. Critically, we render the inversion scheme blind during half a period of the oscillation. We then invert the model (under AR priors on hidden states), with and without large backward lag.