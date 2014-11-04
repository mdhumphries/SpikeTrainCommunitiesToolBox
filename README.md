SpikeTrainCommunitiesToolBox
============================

MATLAB toolbox for spike-train community detection

A set of functions for analysing large-scale recordings of cellular-level neural activity, based on community detection ideas from network theory.

Large-scale recording technology for single-cell activity is now routinely available in variety of methods: silicon probes, multi-electrode arrays, tetrodes, calcium imaging, and voltage-sensitive dye imaging. Having captured the activity of large populations of neurons at single-cell resolution, the next question is: how do I analyse that data?

Key to that analysis is dimension-reduction. One approach to dimension-reduction is to use the fact that neurons tend to fire together in groups - or "ensembles".

We showed how the idea of community detection on arbitrary networks are ideally suited to solve the problem of detecting neural ensembles (Humphries, 2011). 

The purpose of this toolbox is to develop the community-detection algorithms best-suited for the ensemble-detection problem. 

The original code released with Humphries (2011) has been updated to include "consensus" community detection, that dramaticaly improves the reiability of the clustering. For details, see Documentation.

The toolbox includes code to run the core algorithms for community detection, and top-level helper functions that take arrays of spike-train times as input, and do all processing necessary to run the ensemble-detection algorithms. These are developed separately as the core algorithms can run on any suitable similarity matrix, not just those constructed by the top-level functions.

Original paper (supplied in Documentation/ folder): 
Humphries, M. D. (2011) Spike-train communities: finding groups of similar spike trains J Neurosci, 31, 2321-2336

### Undocumented features of the core algorithms:
(1) Suitable for any network: undirected or directed, weighted or unweighted. For example, if using directed correlations between spike-trains, then the algorithm will still cluster into ensembles defined by that directed correlation.

## Deployment

Download the ZIP file (button on right-hand menu)

Open anywhere on your MATLAB path

Top-level script "Consensus_Cluster_Dataset" shows how to call the main function, and plots the results
