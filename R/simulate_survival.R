#' Function to simulate multiple datasets for survival data
#' This function simulates survival data for the integration method.
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
#' @param c0 Control rate parameter of random survival times; larger values = smaller survival times.
#' @param correlated Designates whether the design matrix columns should be correlated (common in genomic data).
#' @param equalcor Only applies if correlated = TRUE. Designates if the correlations the same in all datasets.
#' @return A list containing a length vector of weights corresponding to parameter noise (weights, 
#' NOT the integration weights); datasets each with the design matrix (x), survival times (y), 
#' and censoring vector (censor); and the vector of true parameters beta0 (par)
#' @export

simulate_survival = function(K = 3, p = 5, n = c(50, rep(200, K - 1)), 
                             beta0 = rnorm(p, 0, 1), K0 = 1, sigma.beta = 1, 
                             mu.x = 0, sigma.x = 1, numfrac = 0,
                             sigma.y = 0, c0 = 1, correlated = F, 
                             equalcor = F, trace = F) {
  
  weights = rep(0, K-1)
  weights[sample(K-1, K0)] = 1 
  w = sample(which(weights == 1), numfrac) 
  weights[w] = weights[w] * 0.8 
  
  sigma.beta = c(0, 1 - weights) * sigma.beta 
  mu.x = rep(1, K) * mu.x
  sigma.x = rep(1, K) * sigma.x
  sigma.y = rep(1, K) * sigma.y 
  
  if(equalcor == T){ 
    orthosandwich = matrix(rnorm(p ^ 2, sd = sigma.x), p)
    varcovmat = cov2cor(t(orthosandwich) %*% orthosandwich)
  }
  
  beta = list()
  datasets = list()
  
  for (i in 1:K) { 
    
    beta[[i]] = beta0 + rnorm(p, 0, 1) * sigma.beta[i] 
    
    if(correlated == F){ 
      
      x = matrix(rnorm(n[i] * length(beta[[i]]), mean = mu.x[i], sd = sigma.x[i]), ncol=length(beta[[i]]))
      
    } else if(correlated == T & equalcor == T){ 
      
      x = mvrnorm(n[i], rep(mu.x[i], p), varcovmat)
      
    } else { 
      
      orthosandwich = matrix(rnorm(p ^ 2, me), p)
      varcovmat = cov2cor(t(orthosandwich) %*% orthosandwich)
      x = mvrnorm(n[i], rep(mu.x[i], p), varcovmat)
      
    }
    
    yvar = rnorm(n[i], 0, 1) * sigma.y[i] 
    
    times = rexp(n[i],c0 * exp(as.numeric(x %*% beta[[i]]) + yvar)) 
    time.censor = rexp(n[i], c0 * exp(beta[[i]][1] * x[,1] + yvar)) 
    
    censorv = ifelse(times < time.censor, 1, 0) 
    y <- ifelse(times < time.censor, times, time.censor) 
    ord = order(y, decreasing = F)  
    
    datasets[[i]] = list(x = x[ord,], y = y[ord], censor = censorv[ord])
  }
  
  return(list(weights = weights, datasets = datasets, par = beta0)) 
  
}


