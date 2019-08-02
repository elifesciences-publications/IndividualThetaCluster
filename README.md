# IndividualThetaCluster
 
This code (Demo_SingleThetaCluster.m includes Step1, Step2, Step3)show how cluster individual theta cycle of LFP recorded in hippocampal CA1 pyrmidal layer into four clusters:slow gamma; median gamma; early high gamma and late high gamma through k-means clustering;

The code corresponds the method we proposed in the publication at elife 2019; https://elifesciences.org/articles/44320

Author: Lu Zhang;math2437@hotmail.com; tested with Matlab 2017b and 2019a;

toolbox needed:
Toolbox 1:buzsaki Lab code for processing hc-11 data; not necessary if use your own
data:https://github.com/buzsakilab/buzcode

Toolbox 2:Comunity Clustering;not necessary if only try k-means
clustering;http://netwiki.amath.unc.edu/GenLouvain/GenLouvain


Step 1 includes loading hc-11 data;please download the data first if not use your own data; https://crcns.org/data-sets/hc/hc-11/about-hc-11

Step 2 prepare data samples for clustering; LFP, filtered LFP in theta band,LFP theta phase is needed;

Step 3-1:k-means clustering and 4 theta-gamma asignment;

Step 3-2:community clustering;
Noted that this process is time coumsing when data sample is large.
In this demo, only the first 5000 samples were used; It cost ~1 hour calculation if use all samples in this demo;
