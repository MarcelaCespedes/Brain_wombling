# Brain_wombling
This repository provides users full access to the brain wombling model described in the book chapter 

Cespedes, M.I., McGree, J.M., Drovandi, C.C., Mengersen, K.L., Reid, L.B., Doecke, J.D., Fripp, J., (2020), Chapter 7 titled 'A Bayesian hierarchical approach to jointly model cortical thickness and covariance networks' in Case studies in applied Bayesian data science, editors Mengersen, K., Pudlo, P., Robert, C.P., France: Springer, pages 155-213. 

link to the book is available [here](https://www.springer.com/gp/book/9783030425524).

# This repository contains 
1. The R code to simulate data, accortding to the methods described in the chapter on 35 regions of interest (ROI).
2. R code to run the Brain wombling algorithm on the simulated data (it takes several hours to run the algorithm on a standard PC).
3. Results from an earlier run, and plots from the processed MCMC chains.

As described in the book chapter, the brain wombling algorithm was applied to longitudinal MRI ROI data from the Australian Imaging, Biomarker's and Lifestyle longitudinal study of ageing ([AIBL](https://aibl.csiro.au/research/support/)) data. Refer to link for request to access data from this study. Alternatively, the Alzheimer's Disease Neuroimaging Initiative (ADNI) is another world class ongoing longitudinal study of ageing, whose objectives include making data available to the scientific community without embargo. For more information see Weiner, Michael W., et al. "The Alzheimer’s Disease Neuroimaging Initiative: A review of papers published since its inception." Alzheimer's & Dementia 8.1 (2012): S1-S68. For data aquisition details see this [link](http://adni.loni.usc.edu/).

# Extension of algorithm
The current algorithm allows for patient specific variables to be included primarily as patient specific fixed effects at the highest level of the hierarchical model. An extension of this work is the _dynamic wombling_ algorithm which allows users to incorporate patient level covariance variables at the third level of the hierarchical model, in addition to the main effects at the highest level of the hierarchical model. The dynamic wombling algorithm is described in detail in this [article](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8568).  

For any information/feedback/bugs or comments on this code or the manuscript, please email:
Marcela.Cespedes@csiro.cau
 
