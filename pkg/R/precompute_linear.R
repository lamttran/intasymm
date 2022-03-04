#' Function to precompute linear regression terms for integration speedup
#' 
#' @param datasets A list of K datasets (each dataset having design matrix X and response y).
#' It is assumed the first indexed dataset is the local dataset.
#' @return Precomputed regression terms (X'X and X'y)
#' @export
#' 

precompute_linear = function(datasets, trace = F) {
  
  K = length(datasets)
  precomputed.data = list()
  
  for (i in 1:K) {
    precomputed.data[[i]] = list(xtx = t(datasets[[i]]$x) %*% datasets[[i]]$x, 
                                 xty = t(datasets[[i]]$x) %*% datasets[[i]]$y)
  }
  
  return(precomputed.data)
  
}