#' Function to conduct a likelihood ratio test (LRT), comparing the local dataset to external datasets
#' 
#' @param datasets A list of datasets (each dataset having design matrix X, response y).
#' It is assumed the first indexed dataset is the local dataset.
#' @return A list of p-values corresponding to testing the null hypothesis that one shared set of parameters
#' is sufficient to describe the local dataset and each external dataset
#' @export

test_coefs_linear = function(datasets, trace = F) {
  
  K = length(datasets)
  p.values = rep(1, K - 1)
  for (i in 1:(K - 1)) {
    
    n1 = length(datasets[[1]]$y)
    ni = length(datasets[[i + 1]]$y)
    p = ncol(datasets[[1]]$x)
    
    y1 = c(datasets[[1]]$y, datasets[[i + 1]]$y) #stacked reduced model
    x1 = rbind(datasets[[1]]$x, datasets[[i + 1]]$x)
    
    lm1 = lm(y1 ~ 0 + x1)
    logLik1 = logLik(lm1)
    
    y2 = c(datasets[[1]]$y, datasets[[i + 1]]$y) #block diagonal full model
    x2 = rbind(cbind(datasets[[1]]$x, matrix(0, n1, p)), cbind(matrix(0, ni, p), datasets[[i + 1]]$x))
    lm2 = lm(y2 ~ 0 + x2)
    logLik2 = logLik(lm2)
    
    LR = 2 * (logLik2 - logLik1) #deviance estimator
    p.values[i] = pchisq(LR, df = p, lower.tail = F)
    
  }
  
  return(list(p.values = p.values))
  
}