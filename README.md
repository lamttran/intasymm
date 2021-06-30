# intasymm
Asymmetric Integration of External Datasets to Small Local Data

Before using, please download the required RcppArmadillo functions from https://github.com/lamttran/intasymmRcppArma

# Try it out
You can quickly run some simulated examples by using the following functions:

- Linear models with continuous data: run_simulation_linear() 
- Cox models with survival data: run_simulation_survival()

Due to limitations in how the pec() function interprets objects, run_simulation_survival() is currently hardcoded to only accept p (the number of parameters) equal to 5. 

Numerical output for estimation and prediction error is outputted as the log of the local-only error divided by the post-integration (or lasso or ridge) weighted error. That is, numbers greater than 0 represent improvement in these errors.
