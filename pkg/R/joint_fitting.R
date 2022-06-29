#' A joint fitting function to calculate integration weights for the full optimization and
#' 2-parameter optimization methods as well as their post-integration weighted parameter
#' estimates. 
#' 
#' @import glmnet
#' @param datasets A list of K datasets (each dataset having design matrix x, response y, and censoring vector censor).
#' It is assumed the first indexed dataset is the local dataset.
#' @param family The model to be fit, one of "linear", "logistic", or "cox"
#' @param method One of:
#' \itemize{
#' \item opt - finds the length K-1 vector of integration weights maximizing the cross-validated partial likelihood
#' \item testingopt - using the better of testingpval or testingbinary as a starting point
#' find the 2-vector of p-value cutpoints that maximize the CVPL, turning a K-1 dimensional problem
#' into a 2 dimensional problem. For further details, refer to the loocv_weighted_survival_neg_2param
#' documentation
#' }
#' @param trace if TRUE, additionally prints family and method along with best weights
#' @return A list of integration weights minimizing the CVPL,
#' the CVPL, the post-integration estimates, and the local-only estimates assuming no integration.
#' @export

joint_fitting = function(datasets, family = c("linear", "logistic", "cox"), 
                         method = c("opt", "testingopt"), trace = F){
  K = length(datasets)
  best.weights = rep(0, K-1) #initialize weights as a vector of 0s
  
  if(family == "linear"){
    fit.weighted = fit_weighted_linear 
    loocv.weighted = loocv_weighted_linear
    precompute = precompute_linear 
    test.coefs = test_coefs_linear
    
    precomputed.data = precompute(datasets)
  } else if(family == "logistic"){
    fit.weighted = fit_weighted_logistic 
    loocv.weighted = loocv_weighted_logistic 
    test.coefs = test_coefs_logistic
    
  } else{
    fit.weighted = fit_weighted_survival 
    loocv.weighted = loocv_weighted_survival 
    test.coefs = test_coefs_survival
  }
  
  if(method == "opt"){
    if(family == "linear"){
      opt = optim(par = rep(0.5, K-1), fn = loocv.weighted, datasets = datasets, 
                  precomputed.data = precomputed.data, method = "L-BFGS-B", 
                  lower = rep(0, K-1), upper = rep(1, K-1))
      
      best.weights = opt$par #set weights and LOOCV as output of optim
      loocv = opt$value
    } else if(family == "logistic"){
      opt = optim(par = rep(0.5, K-1), fn = loocv_weighted_logistic_neg, datasets = datasets, 
                  method = "L-BFGS-B", lower = rep(0, K-1), upper = rep(1, K-1))
      best.weights = opt$par
      loocv = opt$value
    } else{
      opt = optim(par = rep(0.5, K-1), fn = loocv_weighted_survival_neg, datasets = datasets, 
                  method = "L-BFGS-B", lower = rep(0, K-1), upper = rep(1, K-1))
      best.weights = opt$par
      loocv = -(opt$value)
    }
  } else if(method == "testingopt"){
    if(family == "linear"){
      testingopt = optim(par = c(0.05, 0.95), fn = loocv_weighted_2param, datasets = datasets, 
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))$par
      
      pvals = unlist(test.coefs(datasets)$p.values)
      
      for(i in 1:(length(best.weights))){#converting p-values to weights
        
        best.weights[i] = ifelse(pvals[i] < testingopt[1], 0,
                                 ifelse(pvals[i] > testingopt[2], 1, (pvals[i] - testingopt[1])/
                                          (testingopt[2] - testingopt[1])))
        
      }
      
      loocv = loocv.weighted(best.weights, datasets)
    } else if(family == "logistic"){
      testingopt = optim(par = c(0.05, 0.95), fn = loocv_weighted_logistic_neg_2param, datasets = datasets, 
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))$par
      
      pvals = unlist(test.coefs(datasets)$p.values)
      
      for(i in 1:(length(best.weights))){#convert p-values to weights
        
        best.weights[i] = ifelse(pvals[i] < testingopt[1], 0,
                                 ifelse(pvals[i] > testingopt[2], 1, (pvals[i] - testingopt[1])/
                                          (testingopt[2] - testingopt[1])))
        
      }
      
      loocv = loocv.weighted(best.weights, datasets, fit.weighted(best.weights, datasets))
    } else{
      testingopt = optim(par = c(0.05, 0.95), fn = loocv_weighted_survival_neg_2param, datasets = datasets, 
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))$par
      
      pvals = unlist(test.coefs(datasets)$p.values)
      
      for(i in 1:(length(best.weights))){#convert p-values to weights
        
        best.weights[i] = ifelse(pvals[i] < testingopt[1], 0,
                                 ifelse(pvals[i] > testingopt[2], 1, (pvals[i] - testingopt[1])/
                                          (testingopt[2] - testingopt[1])))
        
      }
      
      loocv = loocv.weighted(best.weights, datasets, fit.weighted(best.weights, datasets))
    }
  }
  if (trace) {
    print(family)
    print(method)
    print("best.weights")
    print(best.weights)
  }
  
  if(method == "opt" | method == "testingopt"){
    
    return (list(weights = best.weights, par = fit.weighted(best.weights, datasets), 
                 par.localonly = fit.weighted(rep(0, K-1), datasets), loocv = loocv))
    
  }
  
}