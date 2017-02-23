---
title: "The variational Bayes approach"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Typically, the likelihood function contains high-order interaction terms between subsets of the unknown model parameters $$\vartheta$$ (e.g., because of nonlinearities in the model). This implies that the high-dimensional integrals required for Bayesian [parameter estimation](https://en.wikipedia.org/wiki/Estimation_theory) and [model comparison](https://en.wikipedia.org/wiki/Model_selection) cannot be evaluated analytically. Also, it might be computationally very costly to evaluate them using numerical brute force or [Monte-Carlo sampling](https://en.wikipedia.org/wiki/Monte_Carlo_method) schemes. This motivates the use of [variational approaches](https://en.wikipedia.org/wiki/Calculus_of_variations) to approximate Bayesian inference ([Beal, 2003](https://en.wikipedia.org/wiki/Variational_Bayesian_methods)). In brief, **[variational Bayes](https://en.wikipedia.org/wiki/Variational_Bayesian_methods)** (VB) is an iterative scheme that indirectly optimizes an approximation to both the model evidence $$p(y\mid m)$$ and the posterior density $$p(\vartheta\mid y,m)$$. The key trick is to decompose the log model evidence into:

$$\ln p(y\mid m)=F(q)+D_{KL} (q(\vartheta);p(\vartheta\mid y,m))$$


where $$q(\theta)$$  is any density over the model parameters, $$D_{KL}$$ is the [Kullback-Leibler divergence](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence) and the so-called **free energy** $$F(q)$$  is defined as:

$$F(q)=\langle \ln p(y\mid \vartheta,m) \rangle_q-D_{KL}(q(\vartheta);p(\vartheta\mid m))$$

where the expectation $$\langle\rangle_q$$ is taken under $$q$$. One can see that maximizing the functional $$F(q)$$  with respect to $$q$$ indirectly minimizes the Kullback-Leibler divergence between $$q(\vartheta)$$ and the exact posterior $$p(\vartheta\mid y,m)$$. The decomposition of the log evidence is complete in the sense that if $$q(\vartheta)=p(\vartheta\mid y,m)$$, then $$F(q)=\ln p(y\mid m)$$.

The iterative maximization of free energy is done under simplifying assumptions about the functional
form of $$q$$, rendering $$q$$ an approximate posterior density over model parameters and $$F(q)$$ an approximate log model evidence (actually, a lower bound). Typically, one first partitions the model parameters $$\vartheta$$ into distinct subsets and then assumes that $$q$$ factorizes into the product of the ensuing marginal densities. This assumption of "[mean-field](https://en.wikipedia.org/wiki/Mean_field_theory)" separability effectively replaces stochastic dependencies between model variables by deterministic dependencies between the [moments](https://en.wikipedia.org/wiki/Moment_(mathematics)) of their posterior distributions:

$$ \left.
\begin{array}
qq(\vartheta) & = q(\vartheta_1)q(\vartheta_2) \\
\frac{\partial F}{\partial q} & = 0
 \end{array}
\right\}
\implies \ln q(\vartheta_1) = K + \langle \ln\:p(\vartheta\mid m) + \ln(y\mid \vartheta,m)\rangle_{q(\vartheta_2)}$$

where $$K$$ is a normalization constant and we have used a bi-partition of the parameter space ($$\vartheta=\big\{\vartheta_1,\vartheta_2\big\}$$). Critically, the equation's right-hand term can be broken down into a weighted sum of the moments of the distribution $$q(\vartheta_2)$$. This means that one can update the moments of $$q(\vartheta_1)$$ from those of $$q(\vartheta_2)$$, and reciprocally. The equation above can be generalized to any arbitrary mean-field partition and captures the essence of the variational Bayesian approach. The resulting VB algorithm is amenable to analytical treatment (the free energy optimization is made with respect to the moments of the marginal densities), which makes it generic, quick and efficient.

![]({{ site.baseurl }}/images/wiki/VBA/vb1.jpg)

In addition to the mean-field trick, VBA relies upon a further parametric approximation, which essentially consists in summarizing the marginal posterior by their two first-order moments (mean and variance). This local gaussian approximation is known as the **variational-Laplace approach** ([Friston et al. 2007](https://www.ncbi.nlm.nih.gov/pubmed/17055746)).

The statistical properties of the variational-Laplace approach (e.g., in terms of the quality of the ensuing approximation), for the class of models considered in VBA, are described in [Daunizeau et al. (2009)](http://www.sciencedirect.com/science/article/pii/S0167278909002425).

