#' Simulations to assess the relative performance of the integration methods in reducing
#' estimation and prediction error
#'
#' @import pec
#' @import glmnet
#' @import survival
#' @import prodlim
#' @param K The total number of datasets (including the local dataset, which is assumed to be the first dataset)
#' @param p The number of parameters; for now hardcoded to 5 due to requirements of pec() function
#' @param n The sample size of the datasets, should be of length K
#' @param K0 The number of datasets to integrate wrt the local dataset (aka have beta = beta0)
#' @param sigma.beta Additional noise to add to non-integrated or fractionally integrated datasets
#' @param mu.x Mean of the design matrix
#' @param sigma.x Standard deviation of the design matrix
#' @param numfrac Of the K0 to-be-integrated datasets, how many should be fractionally integrated 
#' (i.e. beta has some noise)
#' @param sigma.y Noise to add to linear predictor in calculating exponentially distributed survival times
#' @param m Number of simulation replicates
#' @param methods Methods to run simulations for; refer to joint_fitting_linear documentation for further
#' information. Methods include those in that function as well as "lasso" and "ridge", using glmnet
#' to minimize the leave-one-out error found by shrinking the parameters of the local dataset.
#' @param heterosk Whether the response should be heteroskedastic, dependent on design matrix
#' @param correlated Designates whether the design matrix columns should be correlated (common in genomic data).
#' @param equalcor Only applies if correlated = TRUE. Designates if the correlations the same in all datasets.
#' @param beta0 The true generating parameters, by default a standard normal
#' @param localsize The size of the local dataset (sampled from the first dataset) to be
#' used as training data 
#' @return A matrix of improvements in estimation and prediction error, each column representing a specified
#' method. Values are presented as the log of the local-only error (no integration) divided by the post-
#' integration error; that is positive values represent improvements by performing integration. 
#' @export

run_simulation_linear = function (model = "linear", K = 5, p = 5, n = c(500, rep(200, K - 1)), 
                                  K0 = 2, sigma.beta = 1, mu.x = 0, sigma.x = 1, sigma.y = 0.5, 
                                  m = 50, numfrac = 0, methods = c("opt", "testingopt"), heterosk = F, 
                                  correlated = F, equalcor = F, trace = F, beta0 = rnorm(p, 0, 1),
                                  localsize = 50) {
  
  L2.localonly = NULL #initializing estimation and prediction error matrices
  L2.weighted = NULL
  L2.localonly.pred = NULL
  L2.weighted.pred = NULL
  
  ptm = proc.time()
  
  for (i in 1:m) {
    print(i)
    if (model == "linear") {
      sim = simulate_linear(K = K, p = p, n = n, K0 = K0, sigma.beta = sigma.beta, 
                            mu.x = mu.x, sigma.x = sigma.x, sigma.y = sigma.y, 
                            numfrac = numfrac, correlated = correlated, equalcor = equalcor,
                            trace = trace, beta0 = beta0, heterosk = heterosk, localsize = localsize)
    }
    
    datasets = sim$datasets #datasets, true parameters, and validation dataset
    truepar = sim$par
    validationset = sim$validationset
    
    L2.weighted.method = NULL
    L2.weighted.method.pred = NULL
    fit.pars = NULL
    
    for (method in methods) { #looping over methods
      if (model == "linear") {
        
        if(method == "opt" | method == "testingopt"){ #use joint fitting for full and 2P optimization
          fit = joint_fitting_linear(datasets, fit_weighted_linear, loocv_weighted_linear, precompute_linear, 
                              test.coefs = test_coefs_linear, method = method, trace = trace)
        }
        
        if(method == "lasso"){ #use glmnet to fit lasso
          
          sim.glm = cv.glmnet(datasets[[1]]$x, datasets[[1]]$y, 
                              family = "gaussian", intercept = FALSE, nfolds = length(datasets[[1]]$y),
                              grouped = F)
          
          #choose lambda minimizing mean leave-one-out error, get fitted parameters
          fit = list(weights = rep(0, p), par = as.vector(coef(sim.glm, s = "lambda.min"))[2:(p + 1)],
                     par.localonly = as.vector(lm(datasets[[1]]$y ~ 0 + 
                                                    datasets[[1]]$x)$coefficients))
        }
        
        #estimation error, RMSE on true parameters
        L2.weighted.method = c(L2.weighted.method, sqrt(mean((fit$par - truepar)^2)))
        
        ##prediction error, RMSE on validation dataset
        L2.weighted.method.pred = c(L2.weighted.method.pred, sqrt(mean((validationset$x %*% fit$par - 
                                                                          validationset$y)^2)))
        fit.pars = rbind(fit.pars, fit$par)
      }
    }
    
    
    L2.localonly = c(L2.localonly, sqrt(mean((fit$par.localonly - truepar)^2)))
    L2.localonly.pred = c(L2.localonly.pred, sqrt(mean((validationset$x %*% fit$par.localonly -
                                                          validationset$y)^2)))
    
    L2.weighted = rbind(L2.weighted, L2.weighted.method)
    L2.weighted.pred = rbind(L2.weighted.pred, L2.weighted.method.pred)
    
    if (trace == T) {
      trace.table = rbind(truepar, fit$par.localonly, fit.pars, par.oracle)
      rownames(trace.table) = c("true", "localonly", paste("weighted", methods, sep="."))
      print(trace.table)
    }
  }
  
  colnames(L2.weighted) = paste("w", methods, sep=".")
  colnames(L2.weighted.pred) = paste("w", methods, sep=".")
  
  
  ptm1 = proc.time() - ptm
  
  return(list(PredError = zapsmall(log(L2.localonly.pred/L2.weighted.pred)), 
              EstError = zapsmall(log(L2.localonly/L2.weighted))))
}