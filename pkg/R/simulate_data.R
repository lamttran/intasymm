#' Function to simulate continuous, binary, or survival datasets
#' 
#' @param K The total number of datasets (including the local dataset, which is assumed to be the first dataset)
#' @param p The number of parameters
#' @param n The sample size of the datasets, should be of length K
#' @param beta0 The "true" generating parameters
#' @param c0 For survival data, governs times; greater c0 = shorter times
#' @param K0 The number of datasets to integrate wrt the local dataset (aka have beta = beta0)
#' @param sigma.beta Additional noise to add to non-integrated or fractionally integrated datasets
#' @param mu.x Mean of the design matrix
#' @param sigma.x Standard deviation of the design matrix
#' @param numfrac Of the K0 to-be-integrated datasets, how many should be fractionally integrated (i.e. beta has some noise)
#' @param sigma.y Noise to add to linear predictor
#' @param correlated Designates whether the design matrix columns should be correlated (common in genomic data).
#' @param equalcor Only applies if correlated = TRUE. Designates if the correlations the same in all datasets.
#' @param wt.frac The noise to add to parameters of fractionally weighted external datasets: 
#' noise ~ (1 - wt.frac) * N(0, 1)
#' @param heterosk Heteroskedasticity: makes response error a function of the design matrix
#' @param family The model to generate data for: linear = continuous data, logistic = binary data, cox = survival data
#' @return A list containing a length vector of weights corresponding to parameter noise (weights, 
#' NOT the integration weights); datasets each with the design matrix (x), response (y), 
#' and for family "cox", censoring vector (censor); and the vector of true parameters beta0 (par)
#' @export

simulate_data = function(K = 3, p = 5, n = c(100, rep(200, K - 1)), beta0 = rnorm(p, 0, 1), c0 = 1,
                         K0 = 1, sigma.beta = 1, mu.x = 0, sigma.x = 1, numfrac = 0, sigma.y = 1, 
                         correlated = F, equalcor = F, wt.frac = 0.8, heterosk = F,
                         family = c("linear", "logistic", "cox")){
  weights = rep(0, K-1)
  weights[sample(K-1, K0)] = 1 #randomly choose a dataset to have weight 1
  w = sample(which(weights == 1), numfrac) #numfrac datasets have fractional weight
  weights[w] = weights[w] * wt.frac #these datasets should be partially integrated
  
  sigma.beta = c(0, 1-weights) * sigma.beta
  mu.x = rep(1, K) * mu.x
  sigma.x = rep(1, K) * sigma.x
  sigma.y = rep(1, K) * sigma.y
  
  #because x ~ N(0, 1), covariance and correlation are equal
  #if equalcor == T, the design matrices of all datasets have the same correlation matrix (varcovmat)
  if(equalcor == T){ 
    orthosandwich = matrix(rnorm(p ^ 2), p)
    varcovmat = cov2cor(t(orthosandwich) %*% orthosandwich)
  }
  
  beta = list()
  datasets = list()
  
  for (i in 1:K) {
    beta[[i]] = beta0 + rnorm(p, 0, 1) * sigma.beta[i]
    
    if(correlated == F){ #design matrices are uncorrelated
      x = matrix(rnorm(n[i] * length(beta[[i]])), ncol=length(beta[[i]]))
    } else if(correlated == T & equalcor == T){ #same cor mat for all datasets
      x = mvrnorm(n[i], rep(0, p), varcovmat)
    } else { #different correlation matrix for each dataset
      orthosandwich = matrix(rnorm(p ^ 2), p)
      varcovmat = cov2cor(t(orthosandwich) %*% orthosandwich)
      x = mvrnorm(n[i], rep(0, p), varcovmat)
    }
    
    #Generating the response variable 
    if(family == "linear"){
      if(heterosk == F){
        y = x %*% beta[[i]] + rnorm(n[i], 0, 1) * sigma.y[i]
      } else{
        y = x %*% beta[[i]] + rnorm(n[i], 0, 0.5 + 1/(1 + exp(-rowSums(x)))) * sigma.y[i]
      }
      datasets[[i]] = list(x = x, y = y)
    } else if(family == "logistic"){
      if(heterosk == F){
        z = x %*% beta[[i]] + rnorm(n[i], 0, 1) * sigma.y[i]
      } else{
        z = x %*% beta[[i]] + rnorm(n[i], 0, 0.5 + 1/(1 + exp(-rowSums(x)))) * sigma.y[i]
      }
      y = rbinom(n[i], 1, 1/(1 + exp(-z)))
      datasets[[i]] = list(x = x, y = y)
    } else {
      if(heterosk == F){
        yvar = rnorm(n[i], 0, 1) * sigma.y[i] #variance of exponential parameter to generate times
      } else{
        yvar = rnorm(n[i], 0, 0.5 + 1/(1 + exp(-rowSums(x)))) * sigma.y[i]
      }
      
      times = rexp(n[i],c0 * exp(as.numeric(x %*% beta[[i]]) + yvar)) #death time
      time.censor = rexp(n[i], c0 * exp(beta[[i]][1] * x[,1] + yvar)) #censoring time
      
      censorv = ifelse(times < time.censor, 1, 0) #censoring indicator
      y <- ifelse(times < time.censor, times, time.censor) #observed time
      ord = order(y, decreasing = F)  #order by INCREASING survival time
      
      datasets[[i]] = list(x = x[ord,], y = y[ord], censor = censorv[ord])
    }
  }
  return(list(weights = weights, datasets = datasets, par = beta0))
  
}