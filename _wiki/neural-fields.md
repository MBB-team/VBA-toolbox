---
title: "Neural fields"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

These models are inspired from statistical physics approaches based upon the notion of mean field, i.e. the idea that interactions within micro-scale ensembles of neurons can be captured by summary statistics (i.e., moments of the relevant distribution). They describe the spatio-temporal response of brain networks to experimental manipulations.

At this macroscopic scale, neural states like mean membrane depolarisation can be regarded as a continuum or field, which is a function of space and time. The spatiotemporal dynamics of neural fields essentially depend on how local ensembles influence each other, through connectivity kernels. The latter quantify the amount of anatomical connections as a function of physical distance between any two points on the field. We refer the interested reader to the demonstration script `demo_2DneuralField.m`.

![]({{ site.baseurl }}/images/wiki/tabs/nf1.jpg)
