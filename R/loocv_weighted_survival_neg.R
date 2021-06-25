#' Function to calculate the negative Cox leave-one-out statistic, the cross-validated partial likelihood (CVPL)
#' This function is passed to the base R function optim (which performs minimization) for the 
#' "full optimization" described in the paper, hence it outputs the negative CVPL 
#' (as the CVPL should be maximized).
#' 
#' @importFrom spatstat.utils revcumsum
#' @import intasymmRcppArma
#' @importFrom expm sqrtm
#' @param weights The vector of integration weights, should be equal to the # of external datasets
#' @param datasets A list of datasets (each dataset having design matrix x, response y, and censoring vector censor).
#' It is assumed the first indexed dataset is the local dataset.
#' @return The negative Cox leave-one-out cross-validation statistic for integrating the 
#' external datasets with the given weights (the CVPL). To be passed to optim(). 
#' @export

loocv_weighted_survival_neg <- function(weights, datasets) {
  
  beta = fit_weighted_survival(weights = weights, datasets = datasets) #get fitted betas
  
  theta = as.numeric(exp(datasets[[1]]$x %*% beta)) #again, i-th contribution to hazard
  hazest = cumsum(datasets[[1]]$censor/revcumsum(theta))  #Breslow estimator of baseline hazard H0
  
  if(datasets[[1]]$censor[1] == 0){
    
    hazest = hazest + min(hazest[hazest > 0])/100  #to ensure stability
    
  }
  
  #construct the weight matrix wrt the full likelihood
  Dlocal = diag(theta * hazest)
  
  #1st deriv of full likelihood = X' * delta
  delta = datasets[[1]]$censor - hazest * theta
  xtdx = t(datasets[[1]]$x) %*% Dlocal %*% datasets[[1]]$x
  
  #as before, add in the contributions of the external datasets to the full Hessian
  for (i in 2:length(datasets)) {
    
    theta = as.numeric(exp(datasets[[i]]$x %*% beta))
    hazest = cumsum(datasets[[i]]$censor / revcumsum(theta))
    
    if(datasets[[i]]$censor[1] == 0){
      
      hazest = hazest + min(hazest[hazest > 0])/100  #to ensure stability
      
    }
    
    D = diag(theta * hazest)
    xtdx = xtdx + weights[i-1] * t(datasets[[i]]$x) %*% D %*% datasets[[i]]$x
    
  }
  
  xtdxsoln = solve(xtdx)
  
  #construct V matrix for leave-one-out param calculation
  VmatDIAG = as.numeric(diag(sqrtm(Dlocal) %*% datasets[[1]]$x %*% 
                               xtdxsoln %*% t(datasets[[1]]$x) %*% sqrtm(Dlocal)))
  
  #initialize LOO params
  looParam = beta - t(t(xtdxsoln %*% t(datasets[[1]]$x)) * as.numeric(delta / (1 - VmatDIAG)))
  n = nrow(datasets[[1]]$x) # new row
  
  #cross validated partial likelihood
  thetafull = exp(datasets[[1]]$x %*% looParam) 
  thetaloo = matrix(thetafull[-seq(1, n * n, n + 1)], n - 1, n) #full and loo hazard
  
  censorfull = matrix(rep(datasets[[1]]$censor, n), ncol = n)
  censorloo = matrix((matrix(rep(datasets[[1]]$censor, length(datasets[[1]]$censor)), 
                             ncol = length(datasets[[1]]$censor)))[-seq(1, n * n, n + 1)], n - 1, n)
  
  return(sum(cvpl(theta = thetaloo, censor = censorloo) - cvpl(theta = thetafull, censor = censorfull)))
}