#' Function to simulate multiple datasets for continuous data
#' This function simulates continuous data for integration with linear models
#' 
#' @param K The total number of datasets (including the local dataset, which is assumed to be the first dataset)
#' @param p The number of parameters
#' @param n The sample size of the datasets, should be of length K
#' @param beta0 The "true" generating parameters
#' @param K0 The number of datasets to integrate wrt the local dataset (aka have beta = beta0)
#' @param sigma.beta Additional noise to add to non-integrated or fractionally integrated datasets
#' @param mu.x Mean of the design matrix
#' @param sigma.x Standard deviation of the design matrix
#' @param numfrac Of the K0 to-be-integrated datasets, how many should be fractionally integrated (i.e. beta has some noise)
#' @param sigma.y Noise to add to linear predictor in calculating exponentially distributed survival times
#' @param correlated Designates whether the design matrix columns should be correlated (common in genomic data).
#' @param equalcor Only applies if correlated = TRUE. Designates if the correlations the same in all datasets.
#' @param localsize The size of the local dataset
#' @param val.noerror Should any response error be in the validation dataset?
#' @param heterosk Heteroskedasticity: makes response error a function of the design matrix
#' @return A list containing a length vector of weights corresponding to parameter noise (weights, 
#' NOT the integration weights); datasets each with the design matrix (x), survival times (y), 
#' and censoring vector (censor); and the vector of true parameters beta0 (par)
#' @export

simulate_linear = function(K = 5, p = 5, n = c(500, rep(200, K - 1)), beta0 = rnorm(p, 0, 1), 
                           K0 = 2, sigma.beta = 1, mu.x = 0, sigma.x = 1, sigma.y = 0.5, 
                           numfrac = 0, correlated = F, localsize = 50,
                           equalcor = F, val.noerror = F, heterosk = F, trace = F) {
  
  weights = rep(0, K - 1)
  weights[sample(K - 1, K0)] = 1 #datasets that should be fully integrated
  
  w = sample(which(weights == 1), numfrac)
  weights[w] = weights[w] * 0.8 #datasets that should be partially integrated
  
  sigma.beta = c(0, 1 - weights) * sigma.beta #variance of parameters, design matrix, and response
  mu.x = rep(1, K) * mu.x
  sigma.x = rep(1, K) * sigma.x
  sigma.y = rep(1, K) * sigma.y
  
  if(equalcor == T){ #correlation matrix if design matrix has same cor for all datasets
    orthosandwich = matrix(rnorm(p ^ 2), p)
    varcovmat = cov2cor(t(orthosandwich) %*% orthosandwich)
  }
  
  beta = list()
  datasets = list()
  
  for (i in 1:K) {#looping over datasets
    
    beta[[i]] = beta0 + rnorm(p, 0, 1) * sigma.beta[i]
    
    if(correlated == F){
      
      x = matrix(rnorm(n[i] * length(beta[[i]])), ncol=length(beta[[i]]))
      
    } else if(correlated == T & equalcor == T){ ##X columns are correlated, same cor in local and external datasets
      
      x = mvrnorm(n[i], rep(0, p), varcovmat)
      
    } else { ##X columns are correlated, different cor in local and external datasets
      
      orthosandwich = matrix(rnorm(p ^ 2), p)
      varcovmat = cov2cor(t(orthosandwich) %*% orthosandwich)
      x = mvrnorm(n[i], rep(0, p), varcovmat)
      
    }
    if(heterosk == F){
      y = x %*% beta[[i]] + rnorm(n[i], 0, 1) * sigma.y[i]
    } else{
      y = x %*% beta[[i]] + rnorm(n[i], 0, 0.5 + 1/(1 + exp(-rowSums(x)))) * sigma.y[i]
    }
    datasets[[i]] = list(x = x, y = y)
    
  }
  
  trainrows = sample(1:n[1], localsize, replace = F) #dividing into "local" data and validation dataset
  if(val.noerror == T){
    validationset = list(x = datasets[[1]]$x[-trainrows,], 
                         y = datasets[[1]]$x[-trainrows,] %*% beta0)
  } else if (val.noerror == F){
    validationset = list(x = datasets[[1]]$x[-trainrows,], y = datasets[[1]]$y[-trainrows])
  }
  datasets[[1]] = list(x = datasets[[1]]$x[trainrows,], y = datasets[[1]]$y[trainrows,])
  
  return(list(weights = weights, datasets = datasets, par = beta0, validationset = validationset))
  
}