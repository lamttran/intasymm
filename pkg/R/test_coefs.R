#' Standalone function to compare local and external datasets with a likelihood ratio test (LRT)
#' 
#' @import survival
#' @param datasets A list of datasets (each dataset having design matrix "x", response "y", 
#' and potentially censoring vector "censor").
#' It is assumed the first indexed dataset is the local dataset.
#' @param family The model to perform an LRT for, one of "linear", "logistic", or "cox".
#' @return A list of p-values corresponding to testing the null hypothesis that one shared set of parameters
#' is sufficient to describe the local dataset and each external dataset
#' @export

test_coefs = function(datasets, family = c("linear", "logistic", "cox")){
  K = length(datasets)
  p.values = rep(1, K - 1)
  
  if(family == "linear"){
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
  } else if(family == "logistic"){
    for (i in 1:(K - 1)) {
      n1 = length(datasets[[1]]$y)
      ni = length(datasets[[i + 1]]$y)
      p = ncol(datasets[[1]]$x)
      y1 = c(datasets[[1]]$y, datasets[[i + 1]]$y)
      x1 = rbind(datasets[[1]]$x, datasets[[i + 1]]$x)
      glm1 = glm(y1 ~ 0 + x1, family= "binomial")
      logLik1 = logLik(glm1)
      
      y2 = c(datasets[[1]]$y, datasets[[i + 1]]$y)
      x2 = rbind(cbind(datasets[[1]]$x, matrix(0, n1, p)), cbind(matrix(0, ni, p), datasets[[i + 1]]$x))
      glm2 = glm(y2 ~ 0 + x2, family = "binomial")
      logLik2 = logLik(glm2)
      
      LR = 2 * (logLik2 - logLik1) #deviance estimator
      p.values[i] = pchisq(LR, df = p, lower.tail = F)
    }
  } else{
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
  }
  return(list(p.values = zapsmall(p.values)))
}
