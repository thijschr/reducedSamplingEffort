The effect of reduced sampling effort on recovery of gradients in species composition
=====================================================================================

This repository contains the data and R scripts for the analysis of the effect of reduced
local sampling effort on the recovery of gradients in species composition. The analysis is
based on marine benthic macrofaunal samples sampled at 28 sites in the inner Oslofjord,
each site being represented by four van Veen grab samples each. 

The whole repository can be downloaded to your computer as a zip file (click the Download ZIP on the right-hand side). Unzip the file and make the parent directory your working directory in R. The input data files are located in the data directory, scripts can be found in the R directory. 

Reproduce the analysis by running the analysisSamplingEffort.R file, which will source both the dataPreparation.R and functionsSamplingEffort.R files.

Although not tested on other datasets, the analysis should be able to run on similar datasets having samples nested within sites, where each sample is registered as an individual observation.

The accompanying manuscript has been submitted to a peer-reviewed journal.
