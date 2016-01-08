---
title: "Behavioural DCM"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Behavioural Dynamic Causal Modelling -- or bDCM -- aims at decomposing the brain's transformation of stimuli into behavioural outcomes, in terms of the relative contribution of brain regions and their connections. In brief, bDCM places the brain at the interplay between stimulus and behaviour: behavioural outcomes arise from coordinated activity in (hidden) neural networks, whose dynamics are driven by experimental inputs. Estimating neural parameters that control network connectivity and plasticity effectively performs a neurobiologically-constrained approximation to the brain's input–outcome transform. In other words, neuroimaging data essentially serves to enforce the realism of bDCM's decomposition of input–output relationships. In addition, post-hoc artificial lesions analyses allow us to predict induced behavioural deficits and quantify the importance of network features for funnelling input–output relationships. This is important, because this enables one to bridge the gap with neuropsychological studies of brain-damaged patients. 

![]({{ site.baseurl }}/images/wiki/bdcm/bdcm-schema.png)

For further details, we refer the reader to the original Neuroimage paper: [Rigoux & Daunizeau, 2015](http://www.sciencedirect.com/science/article/pii/S1053811915004231){:target="_blank"}.

The goal of this page is to show you how to implement and analyse your own bDCM using the tools included in the VBA-toolbox. 

# Specifying the classical DCM

The core of bDCM is identical to classical DCM. In the next section, we will briefly review how to prepare your DCM before stepping it up to behavioural predictions. We refer the reader to the [DCM wiki]({{ site.baseurl }}/wiki/dcm) for amore detailed explanations on how to implement a DCM in the VBA-toolbox. There is also plenty of information in the SPM documentation to help you design and preprocess a DCM-friendly experiment.

## Getting the BOLD timeseries

The first step is to extract the BOLD timeseries from the VOI recorded during your fMRI experiment. This can be done easily in SPM:

* Locate the central voxel of your VOI. _eg._ by selecting the peak of activation in the cluster of interest.
* Click on the 'eigenvariate' button. 
* Define the ROI geometry.
* Adjust for effects of interest (or no interest) to clean the data from movement/physio artifacts. 
* Save.

This procedure will provide a timeseries describing the 'typical' activity in the ROI. Repeat for each node.

For each subject, you will have to load all those timeseries and store them in the observation matrix:

```matlab
% load eigenvariates
data_ROI_1 = load('VOI_ROI1_1.mat') ;
data_ROI_2 = load('VOI_ROI2_1.mat') ;
data_ROI_3 = load('VOI_ROI3_1.mat') ;

% store in the observation matrix
y_fmri = [data_ROI_1.Y data_ROI_2.Y data_ROI_3.Y ]' ; 
```
 
## Defining the inputs
 
The inputs to the DCM are in general the same or very similar to the (not convolved) regressors used in the first level analysis. They are stored in a matrix ```u``` that usually has the samenumber of columns as the observation matrix (cf [combining the responses](#combining-the-responses)) or more if you define the inputs at the [micro-time resolution]({{ site.baseurl }}/wiki/Controlling-the-inversion-using-VBA-options/#micro-time-resolution).


## Setting the connectivity

You now have to implement the connectivity affecting your nodes. You need for this to specify 4 set of matrices:

* **Static connectivity**
  The matrix ```A``` defines the invariant connectivity scheme. Origin nodes are the columns, target nodes lines.
  
  ```matlab
  % node 2 projects on node 1
  % node 3 projects on node 2
  % node 1 projects on node 3
  A = [0 1 0 ;   
       0 0 1 ;   
	   1 0 0 ] ; 
  ```
  <br/>

* **Modulatory influences**
  The matrices ```B``` define potential psycho-physiological influences. It should an array of matrices similar coded as A. The different elements of the array are modulated by the respective inputs to the system.
  
  ```matlab
  % the second input modulates the connection from node 3 to node 1
  B{2} = [0 0 0 ;
          0 0 0 ;
		  1 0 0 ] ; 
  ```
  <br/>
  
* **Direct inputs**
  The matrix ```C``` defines which input (columns) enters which node (lines):
  
  ```matlab
  C = [1 0 ;
       0 0 ;
	   1 0 ] ; % first input enter nodes 1 and 3 
  ```
  <br/>

* **Quadratic influences**
  The last set of matrices, the array ```D``` explicits interaction effects:
  
  ```matlab
  % node 2 receive the interaction effect of nodes 3 and 1 
  D{3} = [0 0 0 ;   
          1 0 0 ;   
	      0 0 0 ] ; 
  ```
  <br/>
  
  
If you already performed the DCM analysis in SPM, those matrices correspond respectively to `DCM.a`, `DCM.b`, `DCM.c`, and `DCM.d` of the SPM output structure.

# Extending the DCM for behaviour

Now we have implemented our neural dynamics, we can focus on extending the model to include behavioural predictors.

## Combining the responses

We already have the BOLD observations stored in the ```y_fmri```. The full observation matrix ```y``` should also include the behavioural observations (button responses, skin conductance, etc.) ```y_behaviour```. This might be more tricky than it looks as the multiple sources of observations are usually recorded at different rates. You also need to inform toolbox about respective data distribution. 
We refer the reader to the [mixed observations]({{ site.base_url}}/wiki/Multisources) page which covers this issue.

## Setting the neuro-behavioural mapping

In bDCM, the mapping from neural state patterns to behaviour is based on a quadratic expansion similar to the one describing the functional connectivity between nodes. Below, we describe the set of matrices that can be used to map the DCM nodes activity ```x``` to the behavioural predictors ```r``` that will be fitted (after an eventual transformation by the observation function) to the behavioural responses.

### Direct neural mapping

The most intuitive neuro-bahvioural mapping is a linear predictor that directly maps DCM nodes to behavioural predictors:
![direct neural mapping]({{ site.base_url}}/images/wiki/bdcm/mapping_ha.png){:width="50%"} This is implemented in the matrix ```hA```, where columns representent neural states and lines are the different response predictors:
  
```matlab
hA = [ 0 1 0 ;   % the first response is predicted by node 2
       0 0 1 ] ; % the second response is predicted by node 3
```
  
### Modulated neural mapping

If part of the behaviour is explained by external factors for which you don't have neural correlates in your DCM, you can include an input modulation of the direct mapping:
![modulated neural mapping]({{ site.base_url}}/images/wiki/bdcm/mapping_hb.png){:width="50%"} This is implement in the matrix array `hB`:
  
```matlab
% the first response is predicted by the 3rd node modulated by the 2nd input
hB{2} = [ 0 0 1 ;   
          0 0 0 ] ; 
```

### Direct input mapping
Another solution is to predict the behaviour directly from the inputs, without using the neural computations made within the DCM network:
![cheating mapping]({{ site.base_url}}/images/wiki/bdcm/mapping_hc.png){:width="50%"} In this case, you loose all the benefits of the bDCM approach, as the neural and behavioural observations will be independantly predicted. However, such workaround can provide good control models to benchmark your bDCMs.

This type of "mapping" is implemented in the matrix `hC`:
  
```matlab
% the first response is predicted by a linear combination of the first and second inputs
% the second response is predicted by the last input only
hC = [ 1 1 ;   
       0 1 ] ; 
```

### Quadratic neural mapping

Finally, similar to the quadratic effects in the classic DCM, we can design quadratic predictors. This type of mapping captures gating effects, _ie._ two nodes must be coactivated in order to drive the response.
![quadratic mapping]({{ site.base_url}}/images/wiki/bdcm/mapping_hd.png){:width="50%"} This can be defined in the matrix array `hD`:
  
```matlab
% the 2nd response (line) is jointly predicted by the 1st (array index) and 3rd (column) nodes
hD{1} = [ 0 0 0 ;   
          0 0 1 ] ; 
```

# Defining the complete model

A model is defined by an evolution function, an observation function, and a set of priors (see [here]({{ site.base_url}}/wiki/Structure-of-VBA's-generative-model)). The toolbox already includes generic functions that allow you to implement a bDCM in a couple of lines of code once you defined the connectivity matrices. This section will guide you step by step.

## bDCM evolution function

Once you have defined all the connectivity and mapping matrices, the model can be constructed automatically.  Just call the `prepare_fullDCM` function:

```matlab
%- prepare the bDCM structure
options = prepare_fullDCM( ...
  A, B, C, D,     ... % DCM connectivity
  TR,             ... % sampling timestep
  microDT,        ... % micro_resolution
  homogeneous,    ... % are all nodes similar?
  hA, hB, hC, hD, ... % behavioural mapping	
  sources         ... % observation mixture
);
```

This procedure will include in the options structure all the details of the model. All the dynamics of the bDCM (neural trajectories, baloon model, and response predictors) will be computed by the `@f_DCMwHRFext` function that you must specify as the evolution function of your model.

Note that the `@prepare_fullDCM` function will also set default priors for the bDCM.

## bDCM observation function

Remember we [began this tutorial]({{ site.base_url }}/wiki/behavioural-DCM/#combining-the-responses) by defining a general observation matrix `y` containing both the BOLD timeseries and the behavioural responses. As the toolbox can only handle one observation function to map hidden states to the observation, we need to 1) create one aggregated observation function predicting all the observation types, 2) inform thee toolbox how to split the observations into sub-units.

### Mixing the predictors

For a classical DCM, the toolobox already provides an observation function: `@g_HRF3`. It uses the baloon-model hidden states computed by the evolution function to derive BOLD signal prediction, say `g_fmri`. It will also be the first ingredient of our mixed observation function.

The second element is a transformation for the response predictors. As a matter of fact, the response predictors computed by the `@f_DCMwHRFext` evolution function are continuous variable. If we deal with binary data, we need to map those values onto a probability of response between [0 1]. The simplest way in this case is to use a sigmoid function, like `@g_softmax4decoding`, that will provide say a `g_buttonResp` predictor.

Now if you have a DCM with 3 nodes and 2 button responses, the mixed observation should provide as it's output `g`:

```matlab
% prediction is a 5 elements vector;
g = [g_fmri         ; % 3 lines for the 3 nodes
     g_buttonResp1  ; 
     g_buttonResp2] ;
```

You'll find in the toolbox a default bDCM observation function, `@g_DCMwHRFext`, that does exactly that for you.
It will apply the HRF model on the first observations, depending on the number of nodes in your DCM, and treat all the remaining observations lines as binary responses by applying a sigmoid transform. If you want to implement something different, _eg._ real valued responses, you'll have to adapt the `@g_DCMwHRFext` to your case to provide a coherent prediction to all you obsercations.

### Splitting the observations

The only thing you still have to do is to specify which lines in the observation matrix `y` correspond to the BOLD and to the responses respectively (cf. [mixed observations]({{ site.base_url}}/wiki/Multisources)). For the example given above, this should look like:

```matlab
% specify distribution of observations

% - three BOLD timeseries / nodes (gaussian) 
  sources(1) = struct('out',1:3,'type',0);  
  
% - two binary motor responses (binomial) 
  sources(2) = struct('out',4,'type',1);  % left  hand
  sources(3) = struct('out',5,'type',1);  % right hand
```


# bDCM analysis

We now assume that you are able to construct your bDCM model completely and want to apply it to your data. 

## Simulation and estimation

Behavioural DCM is not different from any other model in the toolbox. You can thus use the generic functions `@simulateNLSS` to simulate artifical data and `@VBA_NLStateSpaceModel` to estimate the model's parameters from your recordings.

You can also look at the demo script `demo_negfeedback` that implement the simple two-nodes bDCM descibed in the original paper ([Rigoux & Daunizeau, 2015](http://www.sciencedirect.com/science/article/pii/S1053811915004231){:target="_blank"}):
![bdcm-demo]({{ site.base_url}}/images/wiki/bdcm/bdcm_example.png)

> **Toy bDCM example from [Rigoux & Daunizeau, 2015](http://www.sciencedirect.com/science/article/pii/S1053811915004231){:target="_blank"}**
The 1st input evokes responses throught the action of the node 1. The 2nd input modulates the influence of the negative feedback loop formed by the node 2.

## Artificial lesions

One of the force of bDCM is its ability to predict the effect of brain lesions or dysconnectivity on the behavioural response. To do so, you simply need to "switch off" the connections you want to test, simulate the reduced model, and analyse the predicted behavioural responses.

In the toolbox, you can automatically test the effect of lesioning respectively all of the nodes of the DCM:

```matlab
% posterior and out come from @VBA_NLStateSpaceModel
results = VBA_bDCM_lesion(posterior, out) ;

% behavioural timeseries predicted when node 2 is lesioned
results.lesion(2).y 

```

![bdcm-lesion]({{ site.base_url}}/images/wiki/bdcm/bdcm_lesion.png){:width="95%"}

> **Example of lesion analysis**
If the 1st node is lesioned (left), then the normal behaviour (plain bars) is completly abolished (hatched bars). If the 2nd node is lesioned (right), we observe an increase in the response rate when the modulatory input is on, showing the loss of the feedback process.

## Susceptibility analysis

Another important measure of the behavioural dependency on the neural architecture is what has be coined 'susceptibility analysis'. It is also an perturbation method that provides a score of how much each connection in the network is important to funnel the effect of a specific input.

Again, this can be done very easily in the toolbox:

```matlab
% posterior and out come from @VBA_NLStateSpaceModel
results = VBA_susceptibility(posterior,out);

% scoring results for each input (line) and connection (column)
results.susceptibility.norm

```


![bdcm-susceptibility]({{ site.base_url}}/images/wiki/bdcm/bdcm_susceptibility.png){:width="85%"}

> **Example of susceptibility analysis**
The first input is mainly dependent on the direct access C to the node 1. The modulatory influence of the second input also requires the two connections forming the feedback loop.
