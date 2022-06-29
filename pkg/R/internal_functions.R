#' Internal fitting and leave-one-out functions for joint fitting and simulation functions 
#' 
#' @import survival
#' @import intasymmRcppArma
#' @importFrom  spatstat.utils revcumsum
#' @importFrom expm sqrtm
#' @param weights numeric vector of weights
#' @param alphabeta 2-vector of p-value cutoffs for integration weights
#' @param datasets A list of K datasets (each dataset having design matrix x, response y, 
#' and potenailly censoring vector censor). It is assumed the first indexed dataset is the local dataset.
#' @param betafitted Fitted vector of parameter estimates, likely from fit_weighted functions
#' @param precomputed.data Precomputed data for PRESS statistic calculation, from precompute_linear
#' @param tol Convergence criteria for logistic and Cox families: stops when the vector 2-norm of parameter estimates
#' differs by less than tol
#' @name internal_functions
NULL
#> NULL

#' @rdname internal_functions
precompute_linear = function(datasets, trace = F) {
K = length(datasets)
precomputed.data = list()
for (i in 1:K) {
  precomputed.data[[i]] = list(xtx = t(datasets[[i]]$x) %*% datasets[[i]]$x, 
                               xty = t(datasets[[i]]$x) %*% datasets[[i]]$y)
}
return(precomputed.data)
}

#' @rdname internal_functions
fit_weighted_linear = function(weights, datasets, precomputed.data = NULL) {
  
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

#' @rdname internal_functions
fit_weighted_logistic <- function(weights, datasets, tol = 1e-6, trace = F) {
  beta_hat <- list()
  p = ncol(datasets[[1]]$x)
  beta_hat[[1]] = matrix(rep(0, p))
  theta = exp(datasets[[1]]$x %*% beta_hat[[1]])
  pred.prob = c(theta/(1 + theta))
  
  xtdx = t(datasets[[1]]$x) %*% diag(pred.prob * (1 - pred.prob)) %*% datasets[[1]]$x
  xtyp = t(datasets[[1]]$x) %*% (datasets[[1]]$y - pred.prob)
  ##above automatically assumes the first dataset has weight 1
  
  for (i in 2:length(datasets)) {
    theta = exp(datasets[[i]]$x %*% beta_hat[[1]])
    pred.prob = c(theta/(1 + theta))
    
    xtdx = xtdx + weights[i-1] * t(datasets[[i]]$x) %*% diag(pred.prob * (1 - pred.prob)) %*% datasets[[i]]$x
    xtyp = xtyp + weights[i-1] * t(datasets[[i]]$x) %*% (datasets[[i]]$y - pred.prob)
  }
  
  beta_hat[[2]] = beta_hat[[1]] + (solve(xtdx) %*% xtyp)
  j = 2
  
  while (norm(beta_hat[[j]] - beta_hat[[j - 1]], type="f") > tol) {
    j = j + 1
    theta = exp(datasets[[1]]$x %*% beta_hat[[j - 1]])
    pred.prob = c(theta/(1 + theta))
    
    xtdx = t(datasets[[1]]$x) %*% diag(pred.prob * (1 - pred.prob)) %*% datasets[[1]]$x
    xtyp = t(datasets[[1]]$x) %*% (datasets[[1]]$y - pred.prob)
    ##above automatically assumes the first dataset has weight 1
    
    for (i in 2:length(datasets)) {
      theta = exp(datasets[[i]]$x %*% beta_hat[[j - 1]])
      pred.prob = c(theta/(1 + theta))
      
      xtdx = xtdx + weights[i-1] * t(datasets[[i]]$x) %*% diag(pred.prob * (1 - pred.prob)) %*% datasets[[i]]$x
      xtyp = xtyp + weights[i-1] * t(datasets[[i]]$x) %*% (datasets[[i]]$y - pred.prob)
    }
    
    beta_hat[[j]] <- beta_hat[[j-1]] + solve(xtdx) %*% xtyp
  }
  
  beta = as.numeric(beta_hat[[j]])
  if(trace == T){
    return(beta_hat)
  } else{
    return(beta)
  }
}

#' @rdname internal_functions
fit_weighted_survival = function(weights, datasets, tol = 1e-6, trace = F){
  
  beta_hat <- list()
  p = ncol(datasets[[1]]$x)
  beta_hat[[1]] = matrix(rep(0, p))
  
  theta = as.numeric(exp(datasets[[1]]$x %*% beta_hat[[1]])) #hazard vector
  loglikIter = sum(datasets[[1]]$censor * log(theta)) - 
    sum(datasets[[1]]$censor * log(revcumsum(theta))) #the initial partial likelihood
  
  P = outer(theta, revcumsum(theta), '/') #revcumsum(theta) is the risk vector
  P[upper.tri(P)] <- 0 #P is the matrix of absolute failure probabilities for each subject at 
  #all observed times
  
  #construct the weight matrix, with elements equal to the 2nd derivatives of the log likelihood
  #wrt the log of the hazard
  W <- -P %*% diag(datasets[[1]]$censor) %*% t(P)
  diag(W) <- diag(P %*% diag(datasets[[1]]$censor) %*% t(1-P))
  
  #xtwx is the Hessian of the partial likelihood
  xtwx = t(datasets[[1]]$x) %*% W %*% datasets[[1]]$x
  
  #xtp is the score of the partial likelihood with respect to beta
  xtp = t(datasets[[1]]$x) %*% (datasets[[1]]$censor - P %*% datasets[[1]]$censor)
  
  #iterate over all the external datasets to construct the overall Hessian
  #and score. Then do a N-R update with b(i) = b(i-1)  + (xtwx)^-1 * xtp
  for(i in 2:length(datasets)){
    
    theta = as.numeric(exp(datasets[[i]]$x %*% beta_hat[[1]]))
    P = outer(theta, revcumsum(theta), '/')
    P[upper.tri(P)] <- 0
    W <- -P %*% diag(datasets[[i]]$censor) %*% t(P)
    diag(W) <- diag(P %*% diag(datasets[[i]]$censor) %*% t(1-P))
    
    xtwx = xtwx + weights[i-1] * t(datasets[[i]]$x) %*% W %*% datasets[[i]]$x
    xtp = xtp + weights[i-1] * t(datasets[[i]]$x) %*% 
      (datasets[[i]]$censor - P %*% datasets[[i]]$censor)
    
  }
  
  beta_hat[[2]] = beta_hat[[1]] + (solve(xtwx) %*% xtp)
  
  #checking if update increases log likelihood 
  #if not, halve step length until likelihood increases
  thetahold = as.numeric(exp(datasets[[1]]$x %*% beta_hat[[2]]))
  loglikehold = sum(datasets[[1]]$censor * log(thetahold)) - 
    sum(datasets[[1]]$censor * log(revcumsum(thetahold))) #step log-likelihood
  
  expo = 0 #counter for step-halving
  betaguess = beta_hat[[2]]
  
  #executed if the step likelihood (loglikehold) is less than the original likelihood (loglikIter)
  while(loglikehold < loglikIter){ #continue halving until step likelihood is greater
    
    expo = expo + 1 
    halving = 0.5 ^ expo
    betaguess = halving * beta_hat[[2]] + (1 - halving) * beta_hat[[1]]
    thetahold = as.numeric(exp(datasets[[1]]$x %*% betaguess))
    loglikehold = sum(datasets[[1]]$censor * log(thetahold)) - 
      sum(datasets[[1]]$censor * log(revcumsum(thetahold)))
    
  }
  
  beta_hat[[2]] = betaguess #set the post step-halving estimates as the new parameters
  
  j = 2
  
  #repeat until the difference between successive parameters is small
  while (norm(beta_hat[[j]] - beta_hat[[j-1]],type="f") > tol) {
    
    j = j + 1
    theta = as.numeric(exp(datasets[[1]]$x %*% beta_hat[[j-1]]))
    loglikIter = sum(datasets[[1]]$censor * log(theta)) - 
      sum(datasets[[1]]$censor * log(revcumsum(theta)))
    
    P = outer(theta, revcumsum(theta), '/')
    P[upper.tri(P)] <- 0
    
    W <- -P %*% diag(datasets[[1]]$censor) %*% t(P)
    diag(W) <- diag(P %*% diag(datasets[[1]]$censor) %*% t(1-P))
    
    xtwx = t(datasets[[1]]$x) %*% W %*% datasets[[1]]$x
    xtp = t(datasets[[1]]$x) %*% (datasets[[1]]$censor - P %*% datasets[[1]]$censor)
    
    for (i in 2:length(datasets)) {
      
      theta = as.numeric(exp(datasets[[i]]$x %*% beta_hat[[j-1]]))
      P = outer(theta, revcumsum(theta), '/')
      P[upper.tri(P)] <- 0
      W <- -P %*% diag(datasets[[i]]$censor) %*% t(P)
      diag(W) <- diag(P %*% diag(datasets[[i]]$censor) %*% t(1-P))
      
      xtwx = xtwx + weights[i-1] * t(datasets[[i]]$x) %*% W %*% datasets[[i]]$x
      xtp = xtp + weights[i-1] * t(datasets[[i]]$x) %*% 
        (datasets[[i]]$censor - P %*% datasets[[i]]$censor)
      
    }
    
    beta_hat[[j]] <- beta_hat[[j-1]] + solve(xtwx) %*% xtp
    
    #same step-halving procedure
    thetahold = as.numeric(exp(datasets[[1]]$x %*% beta_hat[[j]]))
    loglikehold = sum(datasets[[1]]$censor * log(thetahold)) - 
      sum(datasets[[1]]$censor * log(revcumsum(thetahold)))
    
    expo = 0
    betaguess = beta_hat[[j]]
    
    while(loglikehold < loglikIter){ 
      
      expo = expo + 1 
      halving = 0.5^expo
      betaguess = halving * beta_hat[[j]] + (1 - halving) * beta_hat[[j - 1]]
      thetahold = as.numeric(exp(datasets[[1]]$x %*% betaguess))
      loglikehold = sum(datasets[[1]]$censor * log(thetahold)) - 
        sum(datasets[[1]]$censor *log(revcumsum(thetahold)))
      
    }
    
    beta_hat[[j]] = betaguess
    
  }
  
  #return the final update as the fitted beta
  beta = as.numeric(beta_hat[[j]])
  return(beta)
  
}

#' @rdname internal_functions
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

#' @rdname internal_functions
test_coefs_logistic = function(datasets, trace = F) {
  K = length(datasets)
  p.values = rep(1, K - 1)
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
  return(list(p.values = p.values))
}

#' @rdname internal_functions
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

#' @rdname internal_functions
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

#' @rdname internal_functions
loocv_weighted_logistic = function(weights, datasets, betafitted) {
  beta = betafitted
  theta = as.numeric(exp(datasets[[1]]$x %*% beta))
  pred.prob = c(theta/(1 + theta))
  weight.local = diag(pred.prob * (1 - pred.prob))
  xtdx = t(datasets[[1]]$x) %*% weight.local %*% datasets[[1]]$x
  
  for(i in 2:length(datasets)){
    theta = as.numeric(exp(datasets[[i]]$x %*% beta))
    pred.prob = c(theta/(1 + theta))
    xtdx = xtdx + weights[i-1] * t(datasets[[i]]$x) %*% diag(pred.prob * (1 - pred.prob)) %*% datasets[[i]]$x
  }
  
  xtdxsoln = solve(xtdx)
  Dlocalsqrt = sqrtm(weight.local)
  VmatDIAG = diag(Dlocalsqrt %*% datasets[[1]]$x %*% xtdxsoln %*% t(datasets[[1]]$x) %*% Dlocalsqrt)
  theta = as.numeric(exp(datasets[[1]]$x %*% beta))
  pred.prob = c(theta/(1 + theta))
  
  beta.minus = xtdxsoln %*% t(datasets[[1]]$x)
  diff.v = (datasets[[1]]$y - pred.prob)/(1 - VmatDIAG)
  beta.loo = beta - (beta.minus * rep(diff.v, rep(nrow(beta.minus), ncol(beta.minus))))
  rm(beta.minus, diff.v)
  
  theta.loo = diag(exp(datasets[[1]]$x %*% beta.loo))
  pred.prob.loo = c(theta.loo/(1 + theta.loo))
  #we want to minimize MMLCV: mean minus log-likelihood
  MMLCV = -mean(datasets[[1]]$y * log(pred.prob.loo) + (1 - datasets[[1]]$y) * log(1 - pred.prob.loo))
  
  return(MMLCV)
}

#' @rdname internal_functions
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

#' @rdname internal_functions
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

#' @rdname internal_functions
loocv_weighted_logistic_neg = function(weights, datasets) {
  beta = fit_weighted_logistic(weights = weights, datasets = datasets)
  
  theta = as.numeric(exp(datasets[[1]]$x %*% beta))
  pred.prob = c(theta/(1 + theta))
  weight.local = diag(pred.prob * (1 - pred.prob))
  xtdx = t(datasets[[1]]$x) %*% weight.local %*% datasets[[1]]$x
  
  for(i in 2:length(datasets)){
    theta = as.numeric(exp(datasets[[i]]$x %*% beta))
    pred.prob = c(theta/(1 + theta))
    xtdx = xtdx + weights[i-1] * t(datasets[[i]]$x) %*% diag(pred.prob * (1 - pred.prob)) %*% datasets[[i]]$x
  }
  
  xtdxsoln = solve(xtdx)
  Dlocalsqrt = sqrtm(weight.local)
  VmatDIAG = diag(Dlocalsqrt %*% datasets[[1]]$x %*% xtdxsoln %*% t(datasets[[1]]$x) %*% Dlocalsqrt)
  theta = as.numeric(exp(datasets[[1]]$x %*% beta))
  pred.prob = c(theta/(1 + theta))
  
  beta.minus = xtdxsoln %*% t(datasets[[1]]$x)
  diff.v = (datasets[[1]]$y - pred.prob)/(1 - VmatDIAG)
  beta.loo = beta - (beta.minus * rep(diff.v, rep(nrow(beta.minus), ncol(beta.minus))))
  rm(beta.minus, diff.v)
  
  theta.loo = diag(exp(datasets[[1]]$x %*% beta.loo))
  pred.prob.loo = c(theta.loo/(1 + theta.loo))
  #we want to minimize MMLCV: mean minus log-likelihood
  MMLCV = -mean(datasets[[1]]$y * log(pred.prob.loo) + (1 - datasets[[1]]$y) * log(1 - pred.prob.loo))
  
  return(MMLCV)
}

#' @rdname internal_functions
loocv_weighted_logistic_neg_2param = function(alphabeta, datasets) {
  weights = rep(0, length(datasets) - 1)
  pvals = unlist(test_coefs_logistic(datasets)$p.values)
  
  for(i in 1:length(weights)){ #getting weights from p-values
    weights[i] = ifelse(pvals[i] < alphabeta[1], 0,
                        ifelse(pvals[i] > alphabeta[2], 1, 
                               (pvals[i] - alphabeta[1])/ (alphabeta[2] - alphabeta[1])))
  }
  
  beta = fit_weighted_logistic(weights = weights, datasets = datasets)
  
  theta = as.numeric(exp(datasets[[1]]$x %*% beta))
  pred.prob = c(theta/(1 + theta))
  weight.local = diag(pred.prob * (1 - pred.prob))
  xtdx = t(datasets[[1]]$x) %*% weight.local %*% datasets[[1]]$x
  
  for(i in 2:length(datasets)){
    theta = as.numeric(exp(datasets[[i]]$x %*% beta))
    pred.prob = c(theta/(1 + theta))
    xtdx = xtdx + weights[i-1] * t(datasets[[i]]$x) %*% diag(pred.prob * (1 - pred.prob)) %*% datasets[[i]]$x
  }
  
  xtdxsoln = solve(xtdx)
  Dlocalsqrt = sqrtm(weight.local)
  VmatDIAG = diag(Dlocalsqrt %*% datasets[[1]]$x %*% xtdxsoln %*% t(datasets[[1]]$x) %*% Dlocalsqrt)
  theta = as.numeric(exp(datasets[[1]]$x %*% beta))
  pred.prob = c(theta/(1 + theta))
  
  beta.minus = xtdxsoln %*% t(datasets[[1]]$x)
  diff.v = (datasets[[1]]$y - pred.prob)/(1 - VmatDIAG)
  beta.loo = beta - (beta.minus * rep(diff.v, rep(nrow(beta.minus), ncol(beta.minus))))
  rm(beta.minus, diff.v)
  
  theta.loo = diag(exp(datasets[[1]]$x %*% beta.loo))
  pred.prob.loo = c(theta.loo/(1 + theta.loo))
  #we want to minimize MMLCV: mean minus log-likelihood
  MMLCV = -mean(datasets[[1]]$y * log(pred.prob.loo) + (1 - datasets[[1]]$y) * log(1 - pred.prob.loo))
  
  return(MMLCV)
}

#' @rdname internal_functions
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

#' @rdname internal_functions
loocv_weighted_survival_neg_2param = function(alphabeta, datasets) {
  
  weights = rep(0, length(datasets) - 1)
  pvals = unlist(test_coefs_survival(datasets)$p.values)
  
  for(i in 1:length(weights)){ #getting weights from p-values
    weights[i] = ifelse(pvals[i] < alphabeta[1], 0,
                        ifelse(pvals[i] > alphabeta[2], 1, 
                               (pvals[i] - alphabeta[1]) / (alphabeta[2] - alphabeta[1])))
  }
  
  #the same code as above to calculate the cvpl
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