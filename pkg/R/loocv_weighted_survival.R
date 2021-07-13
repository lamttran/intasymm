#' Function to calculate the Cox leave-one-out statistic, the cross-validated partial likelihood (CVPL)
#' 
#' @importFrom spatstat.utils revcumsum
#' @import intasymmRcppArma
#' @importFrom expm sqrtm
#' @param weights The vector of integration weights, should be equal to the # of external datasets
#' @param datasets A list of datasets (each dataset having design matrix x, response y, and censoring vector censor).
#' It is assumed the first indexed dataset is the local dataset.
#' @param betafitted The post-integration weighted parameter estimates, found using fit_weighted_survival()
#' @return The Cox leave-one-out cross-validation statistic for integrating the 
#' external datasets with the given weights (the CVPL). A higher value indicates a better model fit 
#' @export

loocv_weighted_survival <- function(weights, datasets, betafitted) {
  
  beta = betafitted
  
  #use fitted betas to calculate information on leave-one-out parameters
  
  theta = as.numeric(exp(datasets[[1]]$x %*% beta)) ##hazard vector
  hazest = cumsum(datasets[[1]]$censor / revcumsum(theta))  ###Breslow estimator of baseline hazard
  
  if(datasets[[1]]$censor[1] == 0){#small perturbation of estimator to ensure stability
    
    hazest = hazest + min(hazest[hazest > 0]) / 100  
    
  }
  
  #construct the weight matrix wrt the FULL likelihood
  Dlocal = diag(theta * hazest)
  
  #1st deriv of full likelihood = X' * delta
  delta = datasets[[1]]$censor - hazest * theta
  xtdx = t(datasets[[1]]$x) %*% Dlocal %*% datasets[[1]]$x
  
  #as before, add in the contributions of the external datasets, like with the Hessian
  for (i in 2:length(datasets)) {
    
    theta = as.numeric(exp(datasets[[i]]$x %*% beta))
    hazest = cumsum(datasets[[i]]$censor / revcumsum(theta))
    
    if(datasets[[i]]$censor[1] == 0){
      
      hazest = hazest + min(hazest[hazest > 0]) / 100  #to ensure stability when shortest time is censored
      
    }
    
    D = diag(theta * hazest)
    xtdx = xtdx + weights[i-1] * t(datasets[[i]]$x) %*% D %*% datasets[[i]]$x
    
  }
  
  xtdxsoln = solve(xtdx)
  
  #construct V matrix for leave-one-out parameter calculation
  VmatDIAG = as.numeric(diag(sqrtm(Dlocal) %*% datasets[[1]]$x %*% 
                               xtdxsoln %*% t(datasets[[1]]$x) %*% sqrtm(Dlocal)))
  
  #calculate leave-one-out parameters
  looParam = beta - t(t(xtdxsoln %*% t(datasets[[1]]$x)) * as.numeric(delta / (1 - VmatDIAG)))
  n = nrow(datasets[[1]]$x) # new row
  
  #cross validated partial likelihood
  thetafull = exp(datasets[[1]]$x %*% looParam) #full and leave-one-out hazard vectors
  thetaloo = matrix(thetafull[-seq(1, n * n, n + 1)], n - 1, n)
  
  censorfull = matrix(rep(datasets[[1]]$censor, n), ncol = n) #full and LOO censoring vectors
  censorloo = matrix((matrix(rep(datasets[[1]]$censor, length(datasets[[1]]$censor)), 
                             ncol = length(datasets[[1]]$censor)))[-seq(1, n * n, n + 1)], n - 1, n)
  
  #calculate cross-validated partial likelihood. should be maximized
  return(sum(cvpl(theta = thetafull, censor = censorfull) - cvpl(theta = thetaloo, censor = censorloo)))
  
}