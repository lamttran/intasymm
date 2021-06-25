#' Function to perform weighted fitting of a Cox model 
#' 
#' @importFrom  spatstat.utils revcumsum
#' @param weights Integration weights given to datasets
#' @param datasets A list of datasets (each dataset having design matrix x, response y, and censoring vector censor).
#' It is assumed the first indexed dataset is the local dataset.
#' @return A set of fitted post-integration Cox model parameters 
#' @export

fit_weighted_survival = function(weights, datasets, trace = F){
  
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
  while (norm(beta_hat[[j]] - beta_hat[[j-1]],type="f") > 1e-5) {
    
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