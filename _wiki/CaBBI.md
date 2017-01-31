---
title: "Calcium imaging analysis using biophysical models"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


This page describes the demo codes of CaBBI ([Rahmati et al. 2016](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004736), which is a novel method for analyzing (deconvolving) [calcium imaging](https://en.wikipedia.org/wiki/Calcium_imaging) data. More precisely, CaBBI infers neuronal dynamics and parameters from calcium imaging data. For an example of implementation, see the `demo_CaIBB_FHN` and `demo_CaIBB_QGIF` functions.

# Calcium imaging technique
Calcium imaging has been designed to monitor the neuronal activity of both individual neurons and neuronal populations. Cakcium imaging data are fluorescence traces that encode the change in (intracellular) calcium concentration ([Ca2+]). The firing activities are encoded as fluorescence transients, composed of (usually) a rising phase and a decaying phase. In brief, when a neuron fires, [Ca2+] elevates (rising phase) due to the opening of voltage-gated calcium channels (mainly L-type channels), which allows for Ca2+ influx to the cytoplasm. This initial increase is then decayed exponentially by removal mechanisms such as pumping and buffering (decaying phase). Typically, one want sto reconstruct firing activity (e.g., spike trains) from fluorescence traces. However, this is not a trivial problem because data is polluted by noise, it has low temporal resolution (in regular recording configurations), and there is an indirect nonlinear relationship between [Ca2+] kinetics and the underlying membrane potentials.

# CaBBI: A novel method for calcium imaging analysis
In [Rahmati et al. (2016)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004736), a new method called CaBBI (**C**alcium imaging **a**nalysis using **B**iophysical models and **B**ayesian **I**nference) was proposed. This method has several remarkable advantages over previously introduced methods such as the template matching and machine learning methods. The main advantages are as follows. First, CaBBI incorporates biophysical mechanisms of spike and/or burst generation, by using the spiking and bursting neuron models. Using such priors allows for a much more reliable inference about firing events. Second, it can intrinsically be adapted to fluorescence transients with different rise and/or decay times. In particular, it can even reconstruct the firing events from relatively slowly rising fluorescence transients, as with new genetically-encoded calcium indicators. Third, CaBBI not only reconstructs the firing patterns, but also estimates interpretable biophysical states such as membrane potential, [Ca2+] kinetics, and some voltage-gated currents/conductances. For more details, see ([Rahmati et al. 2016](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004736).
CaBBI is based on a biophysical generative models of fluorescence data (see Figure 1) that includes both evolution and observation mappings. The evolution function governing the dynamics of spiking states, calcium channels, and [Ca2+] kinetics. The generated [Ca2+] kinetics are then mapped non-linearly to the fluorescence kinetics through the observation equation. These mappings are used to recover the hidden states given a measured fluoresce trace (through VBA model inversion).

![]({{ site.baseurl }}/images/wiki/CaBBI/Fig1.png)

> **Schena of the CaBBI** Graph illustrating CaBBI's generative model and its inversion. The represented hierarchy in the graph displays how neuronal dynamics relate to fluorescence traces.

CaBBI can be used with different spiking models. In what follows, we demonstrate the well-known Fitzhugh-Nagumo (FHN) model and the so-called "Quadratic-Gaussian Integrate-and-Fire model" (QGIF), which does not require any reset condition (neither for producing the upstroke phase of the spike nor for its re-polarization phase).


## CaBBI demo codes step by step
The demonstration scripts `demo_CaIBB_FHN` and `demo_CaIBB_QGIF` show how to invert the FHN and QGIF generative models, respectively, given a measured fluorescence trace. By default, the demonstration scripts automatically [download the clacium imaging data](https://figshare.com/s/e524c1d214d411e5869c06ec4b8d1f61) that was processed in Rahmati et al. (2016). These *in-vitro* data were recorded from *neonatal* hippocampal tissues with a sampling rate of 22.6 Hz. Note that the ensuing fluorescence transients have relatively very slow rising phases, which render model inversions particularly difficult. This is why recorded data are pre-processed (slow temporal drifts removal) prior to model inversion.

The Figure below is an example of the graphical output of `demo_CaIBB_FHN`:

<!-- insert an image -->
![]({{ site.baseurl }}/images/wiki/CaBBI/Sampel_DemoFHN.png)

> First row: the  pre-processed fluorescence trace, Second row: estimated [Ca2+] kinetics, Third row: estimated membrane potentials, Fourth row: estimated spiking event times (obtained by thresholding the inferred membrane potentials, see "step 10" below). Note: all these variables are automatically saved in `CaBBI`.

Below, we comment the demonstration scripts step by step:

- **Step 1**: The downloaded data contains eight traces: pick the trace you want, e.g.:

  ```matlab
% file name of the fluorescence trace
Fluor_trace_name = 'fluorescence_data7';
```

- **Step 2**: specify the sampling rate (in Hz) of the fluorescence trace.

- **Step 3**: polynomial de-trending method (low frequency drift removal). You can select the degree of the polynomial (either 3 or 4) which will be fitted to the data, e.g.:

  ```matlab
% Polynomial degree (4 or 3)
degree = '4';
```

- **Steps 5 to 7**: Specify VBA inversion options (including priors).

- **Step 8**: Specify the "basal" [Ca2+] value (in nM), e.g.:

  ```matlab
% in [nM], the basal [Ca2+] concentration
Calcium_basal = 50;
```

- **Step 9**: specify the detection threshold that is used to extract the spiking event times from the inferred membrane potentials (cf. the dashed line in the figure above; third row), e.g.:

  ```matlab
% the spike detection threshold
Detection_threshold  = 0;
```

- **Step 10**: To deal with neural noise, an artificial refractory period is enforced after each detected spike event (i.e. all up-crossing fluctuations within the refractory period will be discarded). You can specify this refractory period as follows:

  ```matlab
refracPeriod = 6; % here: 6 msec
if (count > 1) && ( (j-indices(count-1)) > ceil(refracPeriod/dt) )
```

- **Step 11**: save the variables of interest. For example, the estimated decay time-constant of the calcium transients (see also Eq. 17 in [Rahmati et al. 2016](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004736) can be accessed as follows:

  ```matlab
% in [sec]
inferred_Tau_Ca = inF.Tau_Ca0 * exp(posterior.muTheta(1)) * (1/sampling_rate)/dt
```


## You want to adapt CaBBI demos to your data?

Simply specify the directory where your fluorescence trace is saved, as well as your data sampling rate (in Hz). Running the demo scripts will automatically remove the low frequency drifts from your data, and invert the corresponding generative model. The parameters and states estimates will be stored in `CaBBI` (see Step 11).


# References

[Rahmati V, Kirmse K, Marković D, Holthoff K, Kiebel SJ (2016) Inferring Neuronal Dynamics from Calcium Imaging Data Using Biophysical Models and Bayesian Inference. PLoS Comput Biol 12(2): e1004736. doi:10.1371/journal.pcbi.1004736](http://journals.plos.org/ploscompbiol/article?id=info:doi/10.1371/journal.pcbi.1004736){:target="_blank"}


*Author*: [Vahid Rahmati](https://www.researchgate.net/profile/Vahid_Rahmati3){:target="_blank"}, (vahidrahmati92@gmail.com)
