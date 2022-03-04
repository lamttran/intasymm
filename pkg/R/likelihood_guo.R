#' The method of Guo et al., likelihood data fusion, implemented in R for comparison to our method
#' 
#' @param datasets The list of K datasets (each dataset having design matrix x, response y)
#' It is assumed the first indexed dataset is the local dataset.
#' 
#' @return If trace is off, returns parameters estimated by their method. 
#' If trace is on, returns the above, as well as their integraiton weights and 
#' individual coefficients for each dataset
#' @export

likelihood.guo = function(datasets, trace= F){
  K = length(datasets)
  p = ncol(datasets[[1]]$x)
  gammas = rep(1, K)
  weights = rep(1, K)
  x.squared = rep(0, K)
  coeffs.mat = matrix(0, nrow = K, ncol = p)
  
  lm.local = lm(datasets[[1]]$y ~ 0 + datasets[[1]]$x)
  var.local = (summary(lm.local)$sigma) ^ 2
  resids.local = sum((datasets[[1]]$y - datasets[[1]]$x %*% lm.local$coefficients) ^ 2)
  coeffs.mat[1,] = lm.local$coefficients
  x.squared[1] = sum(diag(datasets[[1]]$x %*% t(datasets[[1]]$x)))
  
  for(i in 2:K){
    lm.external = lm(datasets[[i]]$y ~ 0 + datasets[[i]]$x)
    coeffs = lm.external$coefficients
    coeffs.mat[i,] = coeffs
    resids.external = sum((datasets[[1]]$y - datasets[[1]]$x %*% coeffs) ^ 2)
    
    x.squared[i] = sum(diag(datasets[[i]]$x %*% t(datasets[[i]]$x)))
    gammas[i] = exp(1/(2 * var.local) * (resids.local - resids.external))
    
  }
  weights = gammas * x.squared/(sum(gammas * x.squared))
  if(trace == T){
    return(list(gamms = gammas, weights = weights, coeffs.mat = coeffs.mat, 
                par = colSums(weights * coeffs.mat)))
  } else{
    return(par = colSums(weights * coeffs.mat))
  }
}