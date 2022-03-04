# intasymm
Asymmetric Integration of External Datasets to Small Local Data

This package accompanies the paper "A Cross-Validation Statistical Framework for Asymmetric Data Integration" by Lam Tran, Kevin He, Di Wang, and Hui Jiang.

# Starting out
Before using, please install and attach the required RcppArmadillo functions from https://github.com/lamttran/intasymmRcppArma with the following code (the devtools package is required):
- install_github("lamttran/intasymmRcppArma") 
- library(intasymmRcppArma)

Do the same with the intasymm package:
- install_github("lamttran/intasymm", subdir="pkg") 
- library(intasymm)

Please note both packages were both built under R-4.1.0 and users should update to this version of R.

# Try it out with simulated data
After installing and attaching the intasymm package, some simulated examples can be quickly run with the following code:

- Linear models with continuous data: linearsim <- run_simulation_linear()
- Cox models with survival data: coxsim <- run_simulation_survival()

Due to limitations in how the pec() function interprets objects, run_simulation_survival() is currently hardcoded to only accept p (the number of parameters) equal to 5. 

The user can then manipulate the outputs to obtain the following (using linearsim as an example):

- A summary of the relative performance of various integration methods in terms of prediction error or estimation error: 
  - summary(linearsim$PredError) or summary(linearsim$EstError)

Numerical output for estimation and prediction error is outputted as the log of the local-only error divided by the post-integration (or lasso or ridge) weighted error. That is, numbers greater than 0 represent improvement in these errors.

- A boxplot of the percentage reduction in prediction or estimation error (positive values indicating a reduction): 
  - boxplot(100 - 100/exp(linearsim$EstError)) or boxplot(100 - 100/exp(linearsim$PredError))

# Compare with other methods
Two prior methods for weighted data integration are likelihood data fusion, proposed by Guo et al., and the Bayesian power prior, proposed by Ibrahim et al. No R code for likelihood data fusion was available, so we have included their method in this package (likelihood_guo.R). The power prior is implemented in the package NPP. A comparison of our methods with these two other methods is included as the script guo_npp_comparison.R. 

# A real data application
Unfortunately, the SRTR kidney data used in our paper is private and only available after application. Here, we present an alternative real data application using pan-genomic cancer data available at https://github.com/shuanggema/IntePanCancer, specifically their RData file data.Rdata. The process of linking the separate lists of genes, survival times, and censoring vectors via their paper's summary statistics is nontrivial. We therefore include an RData file (entitled cancersets.RData) in the respository for ease of access that was used in our real data analysis.

The RData file contains data from 9 cancer types, each with 7 variables. The first two columns represent the survival time (in months) and the censoring vector (1 = death, 0 = censored) for the patients. The remaining 5 columns are expression data for genes (TMX3, RIBC1, LOC115581, HCCAT5, and LOC644151) chosen by first selecting the 50 most significant genes via univariate Cox models and then choosing the first 5 of these genes to enter a penalized Cox model via glmnet.

We have also included an R script to replicate the real data results. To do so, perform the following steps:
- Attach the intasymmRcppArma and intasymm packages
- Load the cancersets.RData file into the R global environment
- Run the script realdatascript.R

In the interest of time, this script runs through 1 iterate of 2-parameter integration (the paper used 50 iterations of both 2-parameter and full optimization), but the script can be easily extended to fully replicate the results if desired. The raw output used to generate the real data figures are saved as RealDataOutput.xlsx
