#' A joint fitting function to calculate integration weights for the full optimization and
#' 2-parameter optimization methods as well as their post-integration weighted parameter
#' estimates. 
#' 
#' @import glmnet
#' @param datasets A list of K datasets (each dataset having design matrix x, response y, and censoring vector censor).
#' It is assumed the first indexed dataset is the local dataset.
#' @param fit.weighted The fit_weighted_survival function
#' @param loocv.weighted The loocv_weighted_survival function
#' @param test.coefs  The test_coefs_survival function
#' @param method One of:
#' \itemize{
#' \item testingpval - sets the integration weights equal to the likelihood ratio test (LRT) p-values
#' \item testingbinary - if the LRT is significant (p < 0.05), assign integration weight 0; 1 otherwise
#' \item opt - finds the length K-1 vector of integration weights maximizing the cross-validated partial likelihood
#' \item testingopt - using the better of testingpval or testingbinary as a starting point
#' find the 2-vector of p-value cutpoints that maximize the CVPL, turning a K-1 dimensional problem
#' into a 2 dimensional problem. For further details, refer to the loocv_weighted_survival_neg_2param
#' documentation
#' \item lasso - use glmnet to find the lambda minimizing the leave-one-out
#' cross validation error in a lasso-regularized Cox model, using a partial-likelihood loss
#' \item ridge - use glmnet to find the lambda minimizing the leave-one-out
#' cross validation error in a ridge-regularized Cox model, using a partial-likelihood loss
#' }
#' @return For all non lasso/ridge methods, a list of integration weights minimizing the CVPL,
#' the CVPL, the post-integration estimates, and the local-only estimates assuming no integration.
#' For lasso and ridge, all of the above, as well as the penalty value lambda minimizing the
#' leave-one-out error as well as the LOO error value. Note that these methods assume no integration
#' and thus their vector of weights consists of all zeros. 
#' @export

joint_fitting_survival = function(datasets, fit.weighted = fit_weighted_survival, 
                         loocv.weighted = loocv_weighted_survival, 
                         test.coefs = test_coefs_survival, 
                         method = "testingpval", trace = F) {
  
  K = length(datasets)
  best.weights = rep(0, K-1) #initialize weights as a vector of 0s
  
  if (method == "testingpval"){ #identity link between p-values and weights, not used
    
    best.weights = as.numeric(test.coefs(datasets)$p.values)
    loocv = loocv.weighted(best.weights, datasets, fit.weighted(best.weights, datasets))
    
  } else if (method == "opt"){ #full optimization, null hypothesis of no weights
    
    opt = optim(par = rep(0.5, K-1), fn = loocv_weighted_survival_neg, datasets = datasets, 
                method = "L-BFGS-B", lower = rep(0, K-1), upper = rep(1, K-1))
    best.weights = opt$par
    loocv = -(opt$value)
    
  } else if (method == "testingbinary"){ #step cutoff at p = 0.05, not used
    
    best.weights = as.numeric(test.coefs(datasets)$p.values > 0.05)
    loocv = loocv.weighted(best.weights, datasets, fit.weighted(best.weights, datasets))
    
  } else if (method == "testingopt"){ #optimize the 2 param loo, minimizing the loocv value
    ##use the better of binary or pval for the starting point
    
    loocvpval = loocv.weighted(as.numeric(test.coefs(datasets)$p.values), 
                               datasets, fit.weighted(as.numeric(test.coefs(datasets)$p.values), datasets))
    
    loocvbin = loocv.weighted(as.numeric(test.coefs(datasets)$p.values > 0.05), 
                              datasets, fit.weighted(as.numeric(test.coefs(datasets)$p.values > 0.05), 
                                                     datasets))
    if(loocvpval > loocvbin){ #pval is better
      testingopt = optim(par = c(0, 1), fn = loocv_weighted_survival_neg_2param, datasets = datasets, 
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))$par
    } else {
      testingopt = optim(par = c(0.05, 0.05), fn = loocv_weighted_survival_neg_2param, datasets = datasets, 
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))$par
    }
    
    pvals = unlist(test.coefs(datasets)$p.values)
    
    for(i in 1:(length(best.weights))){#convert p-values to weights
      
      best.weights[i] = ifelse(pvals[i] < testingopt[1], 0,
                               ifelse(pvals[i] > testingopt[2], 1, (pvals[i] - testingopt[1])/
                                        (testingopt[2] - testingopt[1])))
      
    }
    
    loocv = loocv.weighted(best.weights, datasets, fit.weighted(best.weights, datasets))
    
  } else if (method == "lasso"){
    #fit a lasso, taking the coefficients minimizing the mean leave-one-out error
    lassofit = cv.glmnet(datasets[[1]]$x, cbind(time = datasets[[1]]$y, status = datasets[[1]]$censor),
                         family = "cox", nfolds = length(datasets[[1]]$y))
    best.weights = rep(0, K-1)
    loocv = loocv.weighted(best.weights, datasets, as.vector(coef(lassofit, s = "lambda.min")))
    
  } else if (method == "ridge"){
    #fit a ridge, taking the coefficients minimizing the mean leave-one-out error
    lassofit = cv.glmnet(datasets[[1]]$x, cbind(time = datasets[[1]]$y, status = datasets[[1]]$censor),
                         family = "cox", nfolds = length(datasets[[1]]$y), alpha = 0)
    best.weights = rep(0, K-1)
    loocv = loocv.weighted(best.weights, datasets, as.vector(coef(lassofit, s = "lambda.min")))
  }
  
  if (trace) {
    print(method)
    print("best.weights")
    print(best.weights)
  }
  
  if(method == "opt" | method == "testingopt" | method == "testingpval" | method == "testingbinary"){
    
    return (list(weights = best.weights, par = fit.weighted(best.weights, datasets), 
                 par.localonly = fit.weighted(rep(0, K-1), datasets), loocv = loocv))
    
  } else if(method == "lasso" | method == "ridge"){ #different because lasso external weights are always 0
    
    return (list(weights = best.weights, par = as.vector(coef(lassofit, s = "lambda.min")), 
                 par.localonly = fit.weighted(rep(0, K-1), datasets), loocv = loocv, 
                 minlambda = lassofit$lambda.min, lambdapos = length(lassofit$lambda) - 
                   which(lassofit$lambda == lassofit$lambda.min)))
    
  }
}