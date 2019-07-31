---
title: "Mixed observations"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

The toolbox allows to develop models of multiple concurrent observations. For example [behavioural-DCM]({{ site.baseurl }}/wiki/behavioural-DCM) can fit both behavioural responses _and_ BOLD timeseries. You may also want to develop a single model explaining both binary button responses _and_ reaction times (e.g., [diffusion-drift models](https://en.wikipedia.org/wiki/Two-alternative_forced_choice#Drift-diffusion_model)). In fact, you can combine virtually as many recordings as you want in a unique generative model where some parameters will impact on multiple observations. In VBA, one can invert such models using the ```sources``` option. The estimation procedure will then try to estimate the parameter values that best explain all data simultaneously.

Please refer to the script ```demo_sources.m``` for an interactive example.

# Creating the observation matrix

To begin with, one needs to form a single observation matrix with all the observations. In what follows, a set of measurements of similar nature (i.e. following a given probability distribution) will be referred to as a "source". For example:

* _button responses_ are a **binary** source. Note: it can also be a multinomial source if more than two alternative button presses are available to subjects within a given trial.
* _(log)-reaction times_ are a **normal** source.
* _BOLD timeseries_ are a **normal** source. If multiple ROIs are recorded, as for a DCM, they can be considered a single multivariate source: VBA will then assume that noise precision is identical accross ROIs.

One forms the multisource observation matrix as follows:

```matlab
y = [ y_source_1 ;
      y_source_2 ;
      ...
      y_source_s ] ;
```

## Case 1: synchronous data

In the simplest case, your observations are always synchronous, i.e. you have one (potentiallay multivariate) observation for each source at each time sample.
Let's say for example you recorded `N` button responses during a learning task, where `N` is the number of trials. You have for each trial ```t``` the observed choice ```y_choice(t)``` and the reaction time ```y_RT(t)```. Then you can construct two horizontal vectors ```y_choice``` and ```y_RT``` having ```N``` values each. The combined observation ```y``` is then directly the 2xN matrix:

```matlab
y = [ y_choice ; y_RT ] ;
```

## Case 2: asynchronous recordings

Different sources can also be recorded with different sampling rates. For example, BOLD timeseries will typically consist of one time sample every 2s or so (1 TR), while button responses will be sparser (with one datapoint per trial, whose onset depends on the ITI, the jitter, and the RT). In such cases, one needs to resample the observations to feed VBa's model inversion with a unique observation matrix, where all datapoints are aligned on the same temporal sampling grid.

First, one chooses a common (re)sampling rate. This sampling rate should be high enough to capture the fastest dynamics. Note that using higher sampling rates will automatically incurr a computational cost (i.e. slower inversion). For behavioural DCM for example, you can resample the data using a sampling period between 10ms and 1s (depending on whether you think trial-by-trial variations in the RTs carry interesting information or not).

First, one resamples the data accordingly. Resampling does not mean interpolation here. Rather, one simply padds the time series with NaNs ("not a number") between datapoints to tell VBA that nothing was recorded there. For example, let's consider an fMRI timeseries that was recorded with a TR of 2s:

```matlab
y_fmri = [ data_mri(1) data_mri(2) ... data_mri(N) ] ;
```

To resample it with a new sampling period of 500ms, we simply insert 3 NaNs between all datapoints:

```matlab
y_fmri_resampled = [data_mri(1) NaN NaN NaN ...
                    data_mri(2) NaN NaN NaN ...
                    ...
                    data_mri(N) NaN NaN NaN] ;
```

A more efficient way to write it is:

```matlab
old_timestep = 2   ; % old TR in s
new_timestep = 0.5 ; % new timestep is 0.5s

% how many new points for an old one
dilution = old_timestep / new_timestep ;

% initialize resampled observation with nans
y_fmri_resampled = nan(1, N*dilution) ;

% insert one mri datapoint every four observation points
y_fmri_resampled(1:dilution:end) = data_mri ;

```

Now, let's say you want to pair the fMRI timeseries with some trial-by-trial button responses (e.g., to perform a behavioural-DCM analysis). Let's assume that you have recorded the onset time ```onset(t)``` of each response ```response(t)```, where the onset times are measured in seconds (beginning at the first fMRI sample). You first have to round those onsets to align them on the same sampling rate as the resampled BOLD. Then, you can use those timings to construct the resampled response vector:

```matlab
% initialize resampled observation with nans
y_resp_resampled = nan(1, numel(y_fmri_resampled)) ;

% find the rounded onsets
onsets_resampled = round(onsets / new_timestep) ;

% set the observed responses at the correct timing
y_resp_resampled(onsets_resampled+1) = response ;

```

It's now time to wrap up. You have one vector for the BOLD timeseries, ```y_fmri_resampled```, and one vector for the behavioural responses, ```y_resp_resampled```, both resampled at the same rate. It is then straightforward to concatenate them to get the mixed observation matrix:

```matlab
y = [y_fmri_resampled ; y_resp_resampled] ;
```

# Defining the mixed observation function

In VBA, [a generative model]({{ site.baseurl }}/wiki/Structure-of-VBA's-generative-model) is specified in terms of both evolution ($$f$$) and observation ($$g$$) functions. Multisource inversion may require to adapt these functions. In particular, one needs to construct an observation function that can provide a prediction for all sources simultanously from the current inputs and hidden states.

Let's assume that we already have an observation function for each observed source, i.e. ```g_source_1(x,u)``` predicts ```y_source_1```, ```g_source_2(x,u)``` predicts ```y_source_2```, and so on. Here, the hidden state ```x``` and the inputs ```u``` must be identical for all sources: a unique (hidden) evolution dynamics underlies all observed time series. Then, the multisource observation function is simply constructed as follows:

```matlab
g(x,u) = [ g_source_1(x,u) ;
           g_source_2(x,u) ;
           ...
           g_source_s(x,u) ] ;
```

For example, behavioural-DCM uses the following multisource observation function:

```matlab
g_bdcm(x,u) = [ g_hrf(x,u)     ;   % predicts BOLD signal
                g_softmax(x,u) ] ; % predicts the binary response
```

# Specifying the observation mixture

We now have a model that can predict multiple sources of observations stored in a unique matrix ```y```. This matrix however contains data following different distributions. In order to derive the data likelihood, we need to specify how each observation should be treated.

To do so, the ```options``` structure must have a ```sources``` field that describes how to split the observation matrix. More precisely, each source ```i``` should be described as follows:

* ```options.sources(i).out``` : the line indices of ```y``` where the corresponding source is stored
* ```options.sources(i).type```: flag for the type of distribution:
	- 0 for (multivariate) gaussian
	- 1 for binomial
	- 2 for for multinomial

In our bDCM example, if we have 3 ROIs/nodes and 1 binary response:

```matlab
% the 3 first lines in the data matrix are normally distributed (BOLD)
options.sources(1).out  = 1:3 ;
options.sources(1).type = 0   ;

% the 4th line in the data matrix is the binary response
options.sources(2).out  = 4 ;
options.sources(2).type = 1   ;

```

# Setting hyperpriors

If the mixed model includes multiple gaussian sources, one should define hyperpriors on noise precision for each of them. This is done by setting the corresponding ```a_sigma``` and ```b_sigma``` fields in the ```options``` structure (see [this page](VBA-model-inversion-in-4-steps) for a complete description of VBA's priors).

To set the hyperpriors for a gaussian source ```s```

```matlab
% Jeffrey's "neutral" hyperprior
options.priors.a_sigma(s) = 1 ;
options.priors.b_sigma(s) = 1 ;
```

Note that the index ```s``` is defined over the gaussian sources and is different from the ```options.sources``` counting. In other words, ```a_sigma``` and ```b_sigma``` should be vectors with a length equal to the total number of gaussian sources only.

When the model is fitted to multiple sources of observation, VBA fuses all information by weighing the different sources according to their respective precision. To reduce potential biases in the ensuing model inversion, one may want to use informed hyperpriors about the sources precisions. This can be done by defining the source-specific hyperprior as a function of the expected range of explained variance:

```matlab
% compute source variance
y_s = y(options.sources(i).out,:);

% set explained variance range
min_ev = 0.05 ; % between 5%
max_ev = 0.15 ; % and 15%

% informed hyperprior
[options.priors.a_sigma(s), options.priors.b_sigma(s)] = ...
	getHyperpriors(y_s,min_ev,max_ev) ;
```

This effectively sets the prior weight of each source in proportion to the amount of expected explained variance.

Multisource inversion can then be performed by calling `VBA_NLStateSpaceModel.m` as usual, but having set the `options` structure as described above.
