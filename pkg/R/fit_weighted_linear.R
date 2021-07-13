#' Function to perform weighted fitting of a linear model 
#' 
#' @param weights Integration weights given to datasets
#' @param datasets A list of datasets (each dataset having design matrix X, response y).
#' It is assumed the first indexed dataset is the local dataset.
#' @param precomputed.data Output from the precompute_linear function if already computed
#' @return A set of fitted post-integration linear model parameters 
#' @export

fit_weighted_linear = function(weights, datasets, precomputed.data = NULL, trace = F) {
  
  if (is.null(precomputed.data)){
    precomputed.data = precompute_linear(datasets)
  }
  
  xtx = precomputed.data[[1]]$xtx
  xty = precomputed.data[[1]]$xty
  
  for (i in 2:length(precomputed.data)) {#weighted sum of x'x, x'y
    xtx = xtx + weights[i-1] * precomputed.data[[i]]$xtx
    xty = xty + weights[i-1] * precomputed.data[[i]]$xty
  }
  
  beta = as.numeric(solve(xtx, xty)) #beta = (x'x)^-1x'y
  return(beta)
}