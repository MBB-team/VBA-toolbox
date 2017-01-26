---
title: "Dynamic Causal Modeling"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Decomposing the relation existing between cognitive functions and their neurobiological "signature" (the spatio-temporal properties of brain activity) requires an understanding of how information is transmitted through brain networks. The ambition here is to ask questions such as: "what is the nature of the information that region A passes on to region B"? This stems from the notion of functional integration, which views function as an emergent property of brain networks. [Dynamic causal modelling](http://www.scholarpedia.org/article/Dynamic_causal_modeling) or DCM was developed specifically to address this question.

DCM embraces a graph-theoretic perspective on brain networks, whereby functionally segregated sources (i.e. brain regions or neuronal populations) correspond to “nodes” and conditional dependencies among the hidden states of each node are mediated by effective connectivity (directed “edges”). DCM generative models are causal in at least two senses:

- DCM describes how experimental manipulations influence the dynamics of hidden (neuronal) states of the system using ordinary differential equations. These evolution equations summarize the biophysical  mechanisms underlying the temporal evolution of states, given a set of unknown evolution parameters that determine both the presence/absence of edges in the graph and how these influence the dynamics of the system’s states.
- DCM maps the system’s hidden states to experimental measures. This observation equation accounts for the main characteristics of the neuroimaging apparatus.

The inversion of such models given neuroimaging data can then be used to identify the structure of brain networks and their specific modulation by the experimental manipulation (i.e. induced plasticity). For example, showing that a given connection is modulated by the saliency of some stimulus demonstrates that this connection conveys the saliency information.

The academic freeware [SPM](http://www.fil.ion.ucl.ac.uk/spm/) proposed the seminal implementations of DCM. The following script allows one to invert an SPM-specified DCM model using VBA, and then eyeball inversion results and diagnostics:

```matlab
[y, u, f_fname, g_fname, dim, options] = dcm2vba(DCM); % SPM-compatible format to VBA-compatible format
[posterior, out] = VBA_NLStateSpaceModel(y, u, f_fname, g_fname, dim, options); % VBA inversion
DCM = vba2dcm(posterior, out, [], TR); % VBA-compatible format to SPM-compatible format
spm_dcm_explore(DCM) % eyeball inversion results
```

where `DCM` is the variable saved in the SPM DCM-file and `TR` is the fMRI repetition time.

![]({{ site.baseurl }}/images/wiki/tabs/dcm1.jpg)

A number of VBA demonstration scripts have been written for DCM:

- `demo_dcm_1region.m`: this script simulates and inverts a 1-node network model. It is useful to assess the statistical identifiability of local self-inhibitory processes and neuro-vascular coupling parameters.
- `demo_dcm4fmri.m`: this script simulates and inverts a 3-nodes network model, with induced network plasticity. Emphasis is put on the influence of neural noise on the system's dynamics.
- `demo_dcm4fmri_distributed.m`: this script simulates and inverts a 3-nodes network model (same as above). Here, the observation function has been augmented with unknown spatial profile of activation, which can be useful to capture non-trivial spatial encoding of experimentally controlled stimuli or observed behaviour.

> **Tip:** By default, a DCM analysis relies upon a deterministic model of neural dynamics. However, VBA can be used to invert *stochastic* DCMs, whereby unpredictable neural perturbations are allowed to interact with responses evoked by the experimental manipulation. This is done by [changing VBA's priors on state noise precision]({{ site.baseurl }}/wiki/VBA-model-inversion-in-4-steps). 

