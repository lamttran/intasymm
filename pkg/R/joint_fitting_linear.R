#' A joint fitting function to calculate integration weights for the full optimization and
#' 2-parameter optimization methods as well as their post-integration weighted parameter
#' estimates. 
#' 
#' @import glmnet
#' @param datasets A list of K datasets (each dataset having design matrix x, response y).
#' It is assumed the first indexed dataset is the local dataset.
#' @param fit.weighted The fit_weighted_linear function
#' @param loocv.weighted The loocv_weighted_linear function
#' @param test.coefs  The test_coefs_linear function
#' @param method One of:
#' \itemize{
#' \item opt - finds the length K-1 vector of integration weights minimizing the PRESS statistic
#' \item testingopt - using the better of testingpval or testingbinary as a starting point
#' find the 2-vector of p-value cutpoints that maximize the CVPL, turning a K-1 dimensional problem
#' into a 2 dimensional problem. For further details, refer to the loocv_weighted_2param
#' documentation
#' }
#' @return A list of integration weights minimizing the PRESS statistic,
#' the PRESS statistic itself, the post-integration estimates, and the 
#' local-only estimates assuming no integration.
#' @export

joint_fitting_linear = function(datasets, fit.weighted = fit_weighted_linear, 
                         loocv.weighted = loocv_weighted_linear, precompute = precompute_linear, 
                         test.coefs = test_coefs_linear, method = c("opt", "testingopt"), trace = F) {
  
  K = length(datasets)
  precomputed.data = precompute(datasets)
  best.weights = rep(0, K-1)
  best.loocv = loocv.weighted(best.weights, datasets, precomputed.data)
  
  if (method == "opt") {#full optimization
    
    opt = optim(par = rep(0.5, K-1), fn = loocv.weighted, datasets = datasets, 
                precomputed.data = precomputed.data, method = "L-BFGS-B", 
                lower = rep(0, K-1), upper = rep(1, K-1))
    
    best.weights = opt$par #set weights and LOOCV as output of optim
    best.loocv = opt$value
    
  } else if (method == "testingopt") {#2 parameter optimization
    
    loocvpval = loocv.weighted(as.numeric(test.coefs(datasets)$p.values), 
                               datasets)
    
    loocvbin = loocv.weighted(as.numeric(test.coefs(datasets)$p.values > 0.05), 
                              datasets)
    
    #set the optimization starting point as the better of identity or binary weighting
    #find the 2 parameters as cutoffs for p-value -> weight function
    
    if(loocvpval < loocvbin){ ##i.e. pval is better = pval has lower loocv error
      
      testingopt = optim(par = c(0, 1), fn = loocv_weighted_2param, datasets = datasets, 
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))$par
    } else {
      
      testingopt = optim(par = c(0.05, 0.05), fn = loocv_weighted_2param, datasets = datasets, 
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))$par
    }
    
    pvals = unlist(test.coefs(datasets)$p.values)
    
    for(i in 1:(length(best.weights))){#converting p-values to weights
      
      best.weights[i] = ifelse(pvals[i] < testingopt[1], 0,
                               ifelse(pvals[i] > testingopt[2], 1, (pvals[i] - testingopt[1])/
                                        (testingopt[2] - testingopt[1])))
      
    }
    
    best.loocv = loocv.weighted(best.weights, datasets)
  }
  
  if (trace) {
    print(method)
    print("best.weights")
    print(best.weights)
  }
  
  return (list(weights = best.weights, par = fit.weighted(best.weights, datasets, precomputed.data), 
               par.localonly = fit.weighted(rep(0, K-1), datasets, precomputed.data), loocv = best.loocv))
}