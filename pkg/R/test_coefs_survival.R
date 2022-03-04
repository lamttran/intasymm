#' Function to conduct a likelihood ratio test (LRT), comparing the local dataset to external datasets
#' 
#' @import survival
#' @param datasets A list of datasets (each dataset having design matrix x, response y, and censoring vector censor).
#' It is assumed the first indexed dataset is the local dataset.
#' @return A list of p-values corresponding to testing the null hypothesis that one shared set of parameters
#' is sufficient to describe the local dataset and each external dataset
#' @export

test_coefs_survival = function(datasets, trace = F) {

K = length(datasets)
p.values = rep(1, K - 1) #vector of LRT p-values

for (i in 1:(K - 1)) { #looping over each external dataset
  
  n1 = length(datasets[[1]]$y)
  ni = length(datasets[[i + 1]]$y)
  p = ncol(datasets[[1]]$x)
  
  #calculating log-likelihood of reduced model
  y1 = c(datasets[[1]]$y, datasets[[i + 1]]$y) #stacked survival times and censoring vectors
  censor1 = c(datasets[[1]]$censor, datasets[[i + 1]]$censor)
  
  x1 = rbind(datasets[[1]]$x, datasets[[i + 1]]$x) #stacked predictors
  cox1 = coxph(Surv(y1, event = censor1) ~ 0 + x1)
  logLik1 = logLik(cox1) #get the log-likelihood of a coxph fit of reduced model
  
  #calculating log-likelihood of a full model
  y2 = c(datasets[[1]]$y, datasets[[i + 1]]$y) #stack the survival times and censoring vectors
  censor2 = c(datasets[[1]]$censor, datasets[[i + 1]]$censor)
  
  x2 = rbind(cbind(datasets[[1]]$x, matrix(0, n1, p)), #make a block diagonal matrix
             cbind(matrix(0, ni, p), datasets[[i + 1]]$x)) #with the local and external datasets
  cox2 = coxph(Surv(y2, event = censor2) ~ 0 + x2)
  logLik2 = logLik(cox2) #get the log-likelihood of a coxph fit of full model
  
  LR = 2 * (logLik2 - logLik1) #deviance estimator
  p.values[i] = pchisq(LR, df = p, lower.tail = F) #perform LRT with p d.f. 
  
}

return(list(p.values = p.values)) #return list of p-values

}