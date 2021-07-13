#' Function to calculate the linear leave-one-out statistic, the PRESS statistic
#' This function is used for full optimization as described in the paper
#' 
#' @param weights The vector of integration weights, should be equal to the # of external datasets
#' @param datasets A list of datasets (each dataset having design matrix X, response y).
#' It is assumed the first indexed dataset is the local dataset.
#' @param precomputed.data Output from the precompute_linear function if already computed
#' @return The linear leave-one-out cross-validation statistic for integrating the 
#' external datasets with the given weights (the PRESS statistic). A lower value indicates a better model fit 
#' @export

loocv_weighted_linear = function(weights, datasets, precomputed.data = NULL) {
  
  if (is.null(precomputed.data)){
    precomputed.data = precompute_linear(datasets)
  }
  
  x = datasets[[1]]$x
  xtx = precomputed.data[[1]]$xtx
  xty = precomputed.data[[1]]$xty
  
  for (i in 2:length(precomputed.data)) {#weighted sum of x'x, x'y
    xtx = xtx + weights[i-1] * precomputed.data[[i]]$xtx
    xty = xty + weights[i-1] * precomputed.data[[i]]$xty
  }
  
  beta = solve(xtx, xty)
  hat = x %*% solve(xtx) %*% t(x)
  yhat = x %*% beta
  loocv = mean(((datasets[[1]]$y - yhat) / (1 - diag(hat)))^2) #calculate PRESS statistic
  
  return(loocv) #Note this statistic should be minimized
}