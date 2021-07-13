#' Function to calculate the linear leave-one-out statistic, the PRESS statistic
#' This function is used for 2 parameter optmization as described in the paper
#' 
#' @param alphabeta The 2-vector of p-value cutpoints; any LRT p-value lower than alphabeta[1]
#' is given an integration weight of 0, any p-value higher than alphabeta[2] is given weight 1
#' @param datasets A list of datasets (each dataset having design matrix X, response y).
#' It is assumed the first indexed dataset is the local dataset.
#' @param precomputed.data Output from the precompute_linear function if already computed
#' @return The linear leave-one-out cross-validation statistic for integrating the 
#' external datasets with the given weights (the PRESS statistic). A lower value indicates a better model fit 
#' @export

loocv_weighted_2param = function(alphabeta, datasets, precomputed.data = NULL){
  
  if (is.null(precomputed.data)){
    precomputed.data = precompute_linear(datasets)
  }
  
  weights = rep(0, length(datasets) - 1)
  pvals = unlist(test_coefs_linear(datasets)$p.values)
  
  for(i in 1:length(weights)){#convert lrt p values to weights
    weights[i] = ifelse(pvals[i] < alphabeta[1], 0,
                        ifelse(pvals[i] > alphabeta[2], 1, 
                               (pvals[i] - alphabeta[1])/ (alphabeta[2] - alphabeta[1])))
  }
  
  x = datasets[[1]]$x
  xtx = precomputed.data[[1]]$xtx
  xty = precomputed.data[[1]]$xty
  
  for (i in 2:length(precomputed.data)) {
    xtx = xtx + weights[i-1] * precomputed.data[[i]]$xtx
    xty = xty + weights[i-1] * precomputed.data[[i]]$xty
  }
  
  beta = solve(xtx, xty)
  hat = x %*% solve(xtx) %*% t(x)
  yhat = x %*% beta
  loocv = mean(((datasets[[1]]$y - yhat) / (1-diag(hat)))^2)
  return(loocv)
  
}