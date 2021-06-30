# intasymm
Asymmetric Integration of External Datasets to Small Local Data

Before using, please download the required RcppArmadillo functions from https://github.com/lamttran/intasymmRcppArma

# Try it out
You can quickly run some simulated examples by using the functions run_simulation_linear (using linear models for continuous data) and run_simulation_survival (using Cox proportional hazards models for survival data). Due to limitations in how the pec() function interprets objects, the survival simulations are currently hardcoded to only accept p (the number of parameters) equal to 5. 

Numerical output for estimation and prediction error is outputted as the log of the local-only error divided by the post-integration (or lasso or ridge) weighted error. That is, numbers greater than 0 represent improvement in these errors.
