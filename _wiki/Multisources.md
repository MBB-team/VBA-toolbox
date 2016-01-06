---
title: "Mixed observations"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

The toolbox allows to develop models of multiple concurrent observations. For example bDCM can fit both behavioural responses _and_ BOLD timeseries. You can also imagine a single model explaining binary button responses _and_ reaction time. In fact, you can combine virtually as many recordings as you want in a unique model where some parameters will affect the prediction of multiple observations. Within the toolbox, this is possible using the ```multisources``` option. The estimation procedure will then try to estimate the parameter values that best explain all data simultaneously.

# Creating the observation matrix

The first step is to create a single observation matrix containing all the observations. We will call 'source' each set of measurements of similar nature, _ie._ following a same probability distribution. For example:

* "_button response_" is a binary source. It can also be a multinomial source if more than two buttons can be choosen from in a given trial.
* "_(log)-reaction time_" is a normally distributed source.
* "_BOLD timeseries_" is a normally distributed source. If multiple ROIs are recorded, as for a DCM, they can be considered a single multivariate source: the recording noise is aproximatively similar accross the regions. 

The general procedure to deal with those multiple observations is first to create the observation matrix:

```matlab
y = [ y_source_1 ;
       y_source_2 ;
       ...
       y_source_n ] ;
```

## Case 1: synchronous data

In the simplest case, your observations are always synchronous, _ie._ you have one (potentiallay multivariate) observation for each source at each timestep. 
Let's say for example you recorded ```N``` button responses in a learning task. You have for each trial ```t``` the selected cue ```y_choice(t)``` and the reaction time ```y_RT(t)```. Then you can construct two horizontal vectors ```y_choice``` and ```y_RT``` having ```N``` values each. The combined observation ```y``` is then directly the 2xN matrix:

```matlab
y = [ y_choice ; y_RT ] ;
```

## Case 2: asynchronous recordings

Different sources can also be recorded with different dynamics. For example, BOLD timeseries will typically consist of one timepoint every 2s or so (1 TR) while button responses will be more sparse with datapoints at times depending on the ITI, the jitter, and the RT. It is then necessary to resample the observations to feed the model with a unique observation matrix where all datapoints are aligned on a same timeline.

The first step is to decide of a common resampling frequency. The new frequency should be high enough to capture the fastest dynamics. Be careful however that the higher the frequency the slower the model estimation. For behavioural DCM for example, you can resample the data using a timestep between 10ms and 1s depending on the importance of the RT. 

The second step is to resample the data. Resampling does not mean interpolation here. It consists more simply into inserting NaNs (not a number) between datapoints to indicate that nothing has been recorded, at the resampling frequency, between those points. For example, let's take an fMRI BOLD timeseries recorded with a TR of 2s:

```matlab
y_fmri = [ data_mri(1) data_mri(2) ... data_mri(N) ] ;
```

If we now want to resample it with a new timestep of 500ms, we should insert 3 NaNs between all datapoints: 

```matlab
y_fmri_resampled = [data_mri(1) NaN NaN NaN ...
                     data_mri(2) NaN NaN NaN ...
                     ... 
                     data_mri(N) NaN NaN NaN] ;
```

A more efficient way to write it is:

```Matlab
old_timestep = 2   ; % old TR in s
 new_timestep = 0.5 ; % new timestep is 0.5s

 % how many new points for an old one
 dilution = old_timestep / new_timestep ; 

 % initialize resampled observation with nans
 y_fmri_resampled = nan(1, N*dilution) ; 
 
 % insert one mri datapoint every four observation points
 y_fmri_resampled(1:dilution:end) = data_mri ; 

```

Now, let's say you want to associate your BOLD signal to some button responses to implement a bDCM. You should have for each response ```response(t)``` an onset time ```onset(t)```, that is the time elapsed between the first volume and the response. You first have to round those onsets to align them on the same sampling frequency as the resampled BOLD. THen, you can use those timings to construct the resampled response vector: 

```matlab
% initialize resampled observation with nans
 y_resp_resampled = nan(1, numel(y_fmri_resampled)) ; 
 
 % find the rounded onsets
 onsets_resampled = round(onsets / new_timestep) ;
 
 % set the observed responses at the correct timing
 y_resp_resampled(onsets_resampled+1) = response ;
 
```

It's now time to wrap up. You have one vector for the BOLD timeseries, ```y_fmri_resampled```, and one vector for the behavioural response, ```y_resp_resampled```, both resampled at the same frequency. It is then straghtforward to concatenate them to get the mixed observation matrix:

```matlab
y = [y_fmri_resampled ; y_resp_resampled] ;
```

# Defining the observation mixture



# Specifying the mixed observation function


# Of the importance of hyperpriors