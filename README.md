# intasymm
Asymmetric Integration of External Datasets to Small Local Data

This package accompanies the paper "A Cross-Validation Statistical Framework for Asymmetric Data Integration" by Lam Tran, Kevin He, and Hui Jiang.

Before using, please install and attach the required RcppArmadillo functions from https://github.com/lamttran/intasymmRcppArma with install_github("lamttran/intasymmRcppArma") and  library(intasymmRcppArma). Please note that package and this one were both built under R-4.1.0 and users should update to this version of R.

# Try it out with simulated data
After installing and attaching the intasymm package, you can quickly run some simulated examples with the following code:

- Linear models with continuous data: linearsim <- run_simulation_linear()
- Cox models with survival data: coxsim <- run_simulation_survival()

Due to limitations in how the pec() function interprets objects, run_simulation_survival() is currently hardcoded to only accept p (the number of parameters) equal to 5. 

The user can then manipulate the outputs to obtain the following (using linearsim as an example):

- A summary of the relative performance of various integration methods in terms of prediction error or estimation error: 
  - summary(linearsim$PredError) or summary(linearsim$EstError)

Numerical output for estimation and prediction error is outputted as the log of the local-only error divided by the post-integration (or lasso or ridge) weighted error. That is, numbers greater than 0 represent improvement in these errors.

- A boxplot of the percentage reduction in prediction or estimation error (positive values indicating a reduction): 
  - boxplot(100 - 100/exp(linearsim$EstError)) or boxplot(100 - 100/exp(linearsim$PredError))


# Replicating the real data analysis in our paper
The cancer real data analysis in our paper used data available at https://github.com/shuanggema/IntePanCancer, specifically their RData file data.Rdata. The process of linking the separate lists of genes, survival times, and censoring vectors via their paper's summary statistics is nontrivial; we therefore include a cleaned version in this package as well (entitled).

