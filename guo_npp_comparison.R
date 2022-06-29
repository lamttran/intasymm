##install tidyr, ggplot2, and NPP if these aren't installed
#install.packages("NPP")
#install.packages("tidyr")
#install.packages("ggplot2")
library(ggplot2)
library(NPP)
library(tidyr)
library(dplyr)

simulate_linear = function(K = 5, p = 5, n = c(500, rep(200, K - 1)), beta0 = rnorm(p, 0, 1), 
                           K0 = 2, sigma.beta = 1, mu.x = 0, sigma.x = 1, sigma.y = 0.5, 
                           numfrac = 0, correlated = F, localsize = 50,
                           equalcor = F, val.noerror = F, heterosk = F, trace = F) {
  
  weights = rep(0, K - 1)
  weights[sample(K - 1, K0)] = 1 #datasets that should be fully integrated
  
  w = sample(which(weights == 1), numfrac)
  weights[w] = weights[w] * 0.8 #datasets that should be partially integrated
  
  sigma.beta = c(0, 1 - weights) * sigma.beta #variance of parameters, design matrix, and response
  mu.x = rep(1, K) * mu.x
  sigma.x = rep(1, K) * sigma.x
  sigma.y = rep(1, K) * sigma.y
  
  if(equalcor == T){ #correlation matrix if design matrix has same cor for all datasets
    orthosandwich = matrix(rnorm(p ^ 2), p)
    varcovmat = cov2cor(t(orthosandwich) %*% orthosandwich)
  }
  
  beta = list()
  datasets = list()
  
  for (i in 1:K) {#looping over datasets
    
    beta[[i]] = beta0 + rnorm(p, 0, 1) * sigma.beta[i]
    
    if(correlated == F){
      
      x = matrix(rnorm(n[i] * length(beta[[i]])), ncol=length(beta[[i]]))
      
    } else if(correlated == T & equalcor == T){ ##X columns are correlated, same cor in local and external datasets
      
      x = mvrnorm(n[i], rep(0, p), varcovmat)
      
    } else { ##X columns are correlated, different cor in local and external datasets
      
      orthosandwich = matrix(rnorm(p ^ 2), p)
      varcovmat = cov2cor(t(orthosandwich) %*% orthosandwich)
      x = mvrnorm(n[i], rep(0, p), varcovmat)
      
    }
    if(heterosk == F){
      y = x %*% beta[[i]] + rnorm(n[i], 0, 1) * sigma.y[i]
    } else{
      y = x %*% beta[[i]] + rnorm(n[i], 0, 0.5 + 1/(1 + exp(-rowSums(x)))) * sigma.y[i]
    }
    datasets[[i]] = list(x = x, y = y)
    
  }
  
  trainrows = sample(1:n[1], localsize, replace = F) #dividing into "local" data and validation dataset
  if(val.noerror == T){
    validationset = list(x = datasets[[1]]$x[-trainrows,], 
                         y = datasets[[1]]$x[-trainrows,] %*% beta0)
  } else if (val.noerror == F){
    validationset = list(x = datasets[[1]]$x[-trainrows,], y = datasets[[1]]$y[-trainrows])
  }
  datasets[[1]] = list(x = datasets[[1]]$x[trainrows,], y = datasets[[1]]$y[trainrows,])
  
  return(list(weights = weights, datasets = datasets, par = beta0, validationset = validationset))
  
}

precompute_linear = function(datasets, trace = F) {
  
  K = length(datasets)
  precomputed.data = list()
  
  for (i in 1:K) {
    precomputed.data[[i]] = list(xtx = t(datasets[[i]]$x) %*% datasets[[i]]$x, 
                                 xty = t(datasets[[i]]$x) %*% datasets[[i]]$y)
  }
  
  return(precomputed.data)
  
}

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

joint_fitting_linear = function(datasets, fit.weighted = fit_weighted_linear, 
                                loocv.weighted = loocv_weighted_linear, precompute = precompute_linear, 
                                test.coefs = test_coefs_linear, method = c("opt", "testingopt"), trace = F) {
  
  K = length(datasets)
  precomputed.data = precompute(datasets)
  best.weights = rep(0, K-1)
  best.loocv = loocv.weighted(best.weights, datasets, precomputed.data)
  
  if (method == "opt") {#full optimization
    
    opt = optim(par = rep(0.5, K-1), fn = loocv.weighted, datasets = datasets, 
                precomputed.data = precomputed.data, method = "L-BFGS-B", 
                lower = rep(0, K-1), upper = rep(1, K-1))
    
    best.weights = opt$par #set weights and LOOCV as output of optim
    best.loocv = opt$value
    
  } else if (method == "testingopt") {#2 parameter optimization
    
    loocvpval = loocv.weighted(as.numeric(test.coefs(datasets)$p.values), 
                               datasets)
    
    loocvbin = loocv.weighted(as.numeric(test.coefs(datasets)$p.values > 0.05), 
                              datasets)
    
    #set the optimization starting point as the better of identity or binary weighting
    #find the 2 parameters as cutoffs for p-value -> weight function
    
    if(loocvpval < loocvbin){ ##i.e. pval is better = pval has lower loocv error
      
      testingopt = optim(par = c(0, 1), fn = loocv_weighted_2param, datasets = datasets, 
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))$par
    } else {
      
      testingopt = optim(par = c(0.05, 0.05), fn = loocv_weighted_2param, datasets = datasets, 
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))$par
    }
    
    pvals = unlist(test.coefs(datasets)$p.values)
    
    for(i in 1:(length(best.weights))){#converting p-values to weights
      
      best.weights[i] = ifelse(pvals[i] < testingopt[1], 0,
                               ifelse(pvals[i] > testingopt[2], 1, (pvals[i] - testingopt[1])/
                                        (testingopt[2] - testingopt[1])))
      
    }
    
    best.loocv = loocv.weighted(best.weights, datasets)
  }
  
  if (trace) {
    print(method)
    print("best.weights")
    print(best.weights)
  }
  
  return (list(weights = best.weights, par = fit.weighted(best.weights, datasets, precomputed.data), 
               par.localonly = fit.weighted(rep(0, K-1), datasets, precomputed.data), loocv = best.loocv))
}

set.seed(2022)

#we used 500 replicates in our paper, for brevity, this example uses 50
m = 50
pred.mat = matrix(0, nrow = 6, ncol = m)
est.mat = matrix(0, nrow = 6, ncol = m)
time.mat = matrix(0, nrow = 6, ncol = m)

for(i in 1:m){
  print(i)
  K = 2
  #These are 3 integration scenarios, with 1 external dataset that should be
  #integrated with varying weights. Uncomment the corresponding lines to compare under
  #different scenarios
  
  #Fully integrate
  sim = simulate_linear(K = K, K0 = 1, n = c(500, rep(200, K - 1)), localsize = 50,
                        sigma.y = c(0.5, rep(0.5, K - 1)), p = 5, beta0 = rnorm(5, 0, 1), numfrac = 0)
  #Partially integrate
  #sim = simulate_linear(K = K, K0 = 1, n = c(500, rep(200, K - 1)), localsize = 50,
  #                      sigma.y = c(0.5, rep(0.5, K - 1)), p = 5, beta0 = rnorm(5, 0, 1), numfrac = 1)
  
  #Don't integrate
  #sim = simulate_linear(K = K, K0 = 0, n = c(500, rep(200, K - 1)), localsize = 50,
  #                      sigma.y = c(0.5, rep(0.5, K - 1)), p = 5, beta0 = rnorm(5, 0, 1), numfrac = 0)
  
  datasets = sim$datasets 
  truepar = sim$par
  validationset = sim$validationset
  
  #full optimization
  ptm = proc.time()
  fit = joint_fitting_linear(datasets, fit_weighted_linear, loocv_weighted_linear, precompute_linear, 
                             test.coefs = test_coefs_linear, method = "opt")
  time.mat[1, i] = (proc.time() - ptm)[3]
  
  #2-par optimization
  ptm = proc.time()
  fit1 = joint_fitting_linear(datasets, fit_weighted_linear, loocv_weighted_linear, precompute_linear, 
                              test.coefs = test_coefs_linear, method = "testingopt")
  time.mat[2, i] = (proc.time() - ptm)[3]
  
  #Guo's likelihood method
  ptm = proc.time()
  fit.guo = likelihood.guo(datasets = datasets, trace = F)
  time.mat[3, i] = (proc.time() - ptm)[3]
  
  #Power prior with NPP
  params_flat = c()
  ptm = proc.time()
  for(j in 2:length(datasets)){
    params_flat = rbind(params_flat, LMNPP_MCMC(y.Cur = datasets[[1]]$y, y.Hist = datasets[[j]]$y, 
                                                x.Cur = datasets[[1]]$x, x.Hist = datasets[[j]]$x, nsample = 1000,
                                                control.mcmc = list(delta.ini = NULL, burnin = 5000, thin = 1))$beta)
  }
  fit.npp.flat = colMeans(params_flat)
  rm(params_flat)
  time.mat[4, i] = (proc.time() - ptm)[3]
  
  #Naive full integration
  ptm = proc.time()
  fit.naive = fit_weighted_linear(weights = rep(1, K - 1), datasets =  datasets)
  time.mat[5, i] = (proc.time() - ptm)[3]
  
  #No integration
  ptm = proc.time()
  fit.zeros = fit_weighted_linear(weights = rep(0, K - 1), datasets =  datasets)
  time.mat[6, i] = (proc.time() - ptm)[3]
  
  #Prediction error matrix
  pred.iter = c(sqrt(mean((validationset$x %*% fit$par - validationset$y)^2)),
                sqrt(mean((validationset$x %*% fit1$par - validationset$y)^2)),
                sqrt(mean((validationset$x %*% fit.guo - validationset$y)^2)),
                sqrt(mean((validationset$x %*% fit.npp.flat[-1] + fit.npp.flat[1] - validationset$y)^2)),
                sqrt(mean((validationset$x %*% fit.naive - validationset$y)^2)),
                sqrt(mean((validationset$x %*% fit.zeros - validationset$y)^2))
  )
  
  pred.mat[,i] = pred.iter
  
  
  #Estimation error matrix
  est.iter = c(sqrt(mean((fit$par - truepar)^2)),
               sqrt(mean((fit1$par - truepar)^2)),
               sqrt(mean((fit.guo - truepar)^2)),
               sqrt(mean((fit.npp.flat[-1] - truepar)^2)),
               sqrt(mean((fit.naive - truepar)^2)),
               sqrt(mean((fit.zeros - truepar)^2))
  )
  
  est.mat[,i] = est.iter
}

est.df = data.frame(t(est.mat))
colnames(est.df) = c("opt", "2par", "guo", "NPP", "naive", "zeros")
est.df = gather(est.df, method, EstError, opt:zeros, factor_key = T)

pred.df = data.frame(t(pred.mat))
colnames(pred.df) = c("opt", "2par", "guo", "NPP", "naive", "zeros")
pred.df = gather(pred.df, method, PredError, opt:zeros, factor_key = T)

time.df = data.frame(t(time.mat))
colnames(time.df) = c("opt", "2par", "guo", "NPP", "naive", "zeros")
time.df = gather(time.df, method, Time, opt:zeros, factor_key = T)

ggplot(est.df, aes(x = method, y = EstError)) + geom_boxplot(fill = "#999999") +
  ggtitle("Estimation Error Comparison") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 

ggplot(pred.df, aes(x = method, y = PredError)) + geom_boxplot(fill = "#999999") + 
  ggtitle("Prediction Error Comparison") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 

ggplot(time.df, aes(x = method, y = Time)) + geom_boxplot(fill = "#999999") + 
  ggtitle("Integration Runtime") + theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = "Runtime (s)") 
