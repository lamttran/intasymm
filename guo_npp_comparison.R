##install tidyr, ggplot2, and NPP if these aren't installed
#install.packages("NPP")
#install.packages("tidyr")
#install.packages("ggplot2")
library(ggplot2)
library(NPP)
library(tidyr)
library(dplyr)

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

ggplot(est.df, aes(x = method, y = EstError)) + geom_boxplot(fill = "green")

ggplot(pred.df, aes(x = method, y = PredError)) + geom_boxplot(fill = "#999999") + 
  ggtitle("External Correct Weight 0") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 

ggplot(time.df, aes(x = method, y = Time)) + geom_boxplot(fill = "#999999") + 
  ggtitle("Integration Runtime") + theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = "Runtime (s)") 
