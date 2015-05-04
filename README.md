## What is the VBA toolbox?

Most models of neurobiological and behavioural data can be broken down into processes that evolve over time and static observation mappings. Given these evolution and observation mappings, the toolbox can be used to simulate data, perform statistical data analysis, optimize the experimental design, etc... In brief, the toolbox provides:

* plug-and-play tools for classical statistical tests
* a library of computational models of behavioural and neurobiological data time series
* quick and efficient probabilistic inference techniques for parameter estimation and model comparison (+ experimental design optimization)
* graphical visualization of results (+ advanced diagnostics of model inversion) 
 
## Requirements
This toolbox runs in Matlab.

## How do I install the toolbox?
1. Get the toolbox
   - [fork](https://github.com/MBB-team/VBA-toolbox/fork) this repo and clone it on your computer ([help!](https://help.github.com/articles/fork-a-repo))
   - download and unzip the ~~[latest stable relase]()~~ or the [current beta version](https://github.com/MBB-team/VBA-toolbox/archive/master.zip)
1. Add the toolbox folder to your Matlab path:
```matlab
cd path_to_/your/VBA_folder
VBA_setup
```
1. Enjoy!
  

## Want more details?
Please visit the [wiki pages](http://mbb-team.github.io/VBA-toolbox/wiki/) for tutorials, demos and advanced features descriptions.

## How can I participate?

This is a collaborative project, in that users can contribute to the toolbox either through [feedback](https://github.com/MBB-team/VBA-toolbox/issues), or directly by including new models for neurobiological and behavioural data. Critically, the toolbox has been developped in the aim of facilitating the inclusion of new models (without having to care about their statistical treatment).
Please send us your code or a [pull request](https://github.com/MBB-team/VBA-toolbox/pulls) so we can include your models in the next release!
