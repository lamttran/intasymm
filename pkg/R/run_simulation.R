#' Simulations to assess the relative performance of the integration methods in reducing
#' estimation and prediction error
#'
#' @import pec
#' @import glmnet
#' @import survival
#' @import prodlim
#' @import pROC
#' @param m The number of replicates 
#' @param family The model to be fit, one of "linear", "logistic", or "cox".
#' @param K The total number of datasets (including the local dataset, which is assumed to be the first dataset)
#' @param p The number of parameters
#' @param beta0 The "true" generating parameters, a vector of length p 
#' @param n The sample size of the datasets, should be of length K
#' @param K0 The number of datasets to integrate wrt the local dataset (aka have beta = beta0)
#' @param sigma.beta Additional noise to add to non-integrated or fractionally integrated datasets
#' @param mu.x Mean of the design matrix
#' @param sigma.x Standard deviation of the design matrix
#' @param sigma.y Noise to add to linear predictor 
#' @param numfrac Of the K0 to-be-integrated datasets, how many should be fractionally integrated 
#' (i.e. beta has some noise)
#' @param wt.frac The noise to add to parameters of fractionally weighted external datasets: 
#' noise ~ (1 - wt.frac) * N(0, 1)
#' @param methods Method(s) to run simulations for; refer to joint_fitting documentation for further
#' information
#' @param heterosk Heteroskedasticity: makes response error a function of the design matrix
#' @param correlated Designates whether the design matrix columns should be correlated (common in genomic data).
#' @param equalcor Only applies if correlated = TRUE. Designates if the correlations the same in all datasets.
#' @param localsize The size of the local dataset (sampled from the first dataset) to be
#' used as training data 
#' @param trace If true, returns additional family and method arguments during joint fitting 
#' @return A matrix of improvements in estimation and prediction error, each column representing a specified
#' method. Values are presented as the log of the local-only error (no integration) divided by the post-
#' integration error; that is positive values represent improvements by performing integration. 
#' @export

run_simulation = function(m = 20, family = c("linear", "logistic", "cox"), K = 5, 
                          p = 5, beta0 = rnorm(p, 0, 1), n = c(500, rep(200, K - 1)), 
                          K0 = 2, sigma.beta = 1, mu.x = 0, sigma.x = 1, 
                          sigma.y = 0, numfrac = 0, wt.frac = 0.8, 
                          methods = c("opt", "testingopt"), heterosk = F, 
                          correlated = F, equalcor = F, localsize = 50, trace = T){
  
  if(family == "linear"){
    PredErrorLocal = NULL #intialize estimation and prediction error matrices
    PredErrorWeighted = NULL
    EstErrorLocal = NULL
    EstErrorWeighted = NULL
    
    for(i in 1:m){
      print(i)
      
      PredErrorweighted.method = NULL #that iteration's prediction and estimation error matrices
      EstErrorweighted.method = NULL
      
      sim = simulate_data(K = K, p = p, n = n, K0 = K0, sigma.beta = sigma.beta, family = "linear",
                          mu.x = mu.x, sigma.x = sigma.x, sigma.y = sigma.y, wt.frac = wt.frac,
                          numfrac = numfrac, correlated = correlated, equalcor = equalcor,
                          beta0 = beta0, heterosk = heterosk)
      
      datasets = sim$datasets #datasets, true parameters
      truepar = sim$par
      trainrows = c(sort(sample(1:nrow(sim$datasets[[1]]$x), size = localsize, replace = FALSE)))
      
      #the validation and local datasets
      local.dataframe.test = list(x = sim$datasets[[1]]$x[-trainrows,], y = sim$datasets[[1]]$y[-trainrows])
      
      local.dataframe.train = list(x = sim$datasets[[1]]$x[trainrows,], y = sim$datasets[[1]]$y[trainrows])
      fit.local = lm(y ~ 0 + x, data = local.dataframe.train)
      sim$datasets[[1]] <- local.dataframe.train
      
      
      fit.pars = NULL
      for(method in methods){
        fit = joint_fitting(datasets, family = "linear", method = method, trace = trace)
        
        EstErrorweighted.method = c(EstErrorweighted.method, sqrt(mean((fit$par - truepar)^2)))
        
        ##prediction error, RMSE on validation dataset
        PredErrorweighted.method = c(PredErrorweighted.method, 
                                     sqrt(mean((local.dataframe.test$x %*% fit$par - 
                                                  local.dataframe.test$y)^2)))
        fit.pars = rbind(fit.pars, fit$par)
      }
      
      EstErrorLocal = c(EstErrorLocal, sqrt(mean((fit$par.localonly - truepar)^2)))
      PredErrorLocal = c(PredErrorLocal, sqrt(mean((local.dataframe.test$x %*% fit$par.localonly -
                                                      local.dataframe.test$y)^2)))
      
      EstErrorWeighted = rbind(EstErrorWeighted, EstErrorweighted.method)
      PredErrorWeighted = rbind(PredErrorWeighted, PredErrorweighted.method)
    }
    
    PredErrorRatio <- log(PredErrorLocal/PredErrorWeighted) #matrices of log pred and est error
    EstErrorRatio <- log(EstErrorLocal/EstErrorWeighted)
    
    colnames(PredErrorRatio) = paste("w", methods, sep=".")
    colnames(EstErrorRatio) = paste("w", methods, sep=".") 
    
    return(list(PredError = zapsmall(PredErrorRatio), #return list of error improvements by method
                EstError = zapsmall(EstErrorRatio)))
    
  } else if(family == "logistic"){
    PredErrorLocal = NULL #intialize estimation and prediction error matrices
    PredErrorWeighted = NULL
    EstErrorLocal = NULL
    EstErrorWeighted = NULL
    
    for(i in 1:m){
      print(i)
      
      PredErrorweighted.method = NULL #that iteration's prediction and estimation error matrices
      EstErrorweighted.method = NULL
      
      sim = simulate_data(K = K, p = p, K0 = K0, n = n, numfrac = numfrac, 
                          sigma.beta = sigma.beta, family = "logistic", sigma.y = sigma.y,
                          wt.frac = wt.frac, correlated = correlated, equalcor = equalcor, 
                          beta0 = beta0, heterosk = heterosk)
      
      trainrows = c(sort(sample(1:nrow(sim$datasets[[1]]$x), size = localsize, replace = FALSE)))
      
      #the validation and local datasets
      local.dataframe.test = list(x = sim$datasets[[1]]$x[-trainrows,], y = sim$datasets[[1]]$y[-trainrows])
      
      local.dataframe.train = list(x = sim$datasets[[1]]$x[trainrows,], y = sim$datasets[[1]]$y[trainrows])
      fit.local = glm(y ~ 0 + x, family = "binomial", data = local.dataframe.train)
      sim$datasets[[1]] <- local.dataframe.train
      
      PredErrorLocal = c(PredErrorLocal, 1 - auc(local.dataframe.test$y, 
                                                 predict(fit.local, local.dataframe.test, type = "response")))
      
      fit.pars = NULL
      for(method in methods){
        fit.method = fit.local
        fit = joint_fitting(datasets = sim$datasets, family = "logistic", method = method, trace = trace) 
        fit.pars = rbind(fit.pars, fit$par)
        
        fit.method$coefficients = fit$par
        
        PredErrorweighted.method = c(PredErrorweighted.method, 1 - 
                                       auc(local.dataframe.test$y, 
                                           predict(fit.method, local.dataframe.test, type = "response")))
        EstErrorweighted.method = c(EstErrorweighted.method, sum((fit$par - sim$par) ^ 2))
      }
      PredErrorWeighted = rbind(PredErrorWeighted, PredErrorweighted.method) #attach est and pred error
      EstErrorLocal = c(EstErrorLocal, sum((fit$par.localonly - sim$par) ^ 2))
      EstErrorWeighted = rbind(EstErrorWeighted, EstErrorweighted.method)
    }
    
    PredErrorRatio <- log(PredErrorLocal/PredErrorWeighted) #matrices of log pred and est error
    EstErrorRatio <- log(EstErrorLocal/EstErrorWeighted)
    
    colnames(PredErrorRatio) = paste("w", methods, sep=".")
    colnames(EstErrorRatio) = paste("w", methods, sep=".") 
    
    return(list(PredError = zapsmall(PredErrorRatio), #return list of error improvements by method
                EstError = zapsmall(EstErrorRatio)))
    
  } else{
    PredErrorLocal = NULL #intialize estimation and prediction error matrices
    PredErrorWeighted = NULL
    EstErrorLocal = NULL
    EstErrorWeighted = NULL
    surv.formula = as.formula(paste("Surv(y, censor) ~ ", paste(paste0("x.", 1:p), collapse= "+")))
    
    for(i in 1:m){
      print(i)
      PredErrorweighted.method = NULL #that iteration's prediction and estimation error matrices
      EstErrorweighted.method = NULL
      
      sim = simulate_data(K = K, p = p, K0 = K0, n = n, numfrac = numfrac, 
                          sigma.beta = sigma.beta, family = "cox", sigma.y = sigma.y,
                          wt.frac = wt.frac, correlated = correlated, equalcor = equalcor, 
                          beta0 = beta0, heterosk = heterosk)
      
      trainrows = c(sort(sample(1:(nrow(sim$datasets[[1]]$x) - 1), size = localsize, 
                                replace = FALSE)), nrow(sim$datasets[[1]]$x)) 
      
      #the local dataset
      local.dataframe.train = data.frame(cbind(sim$datasets[[1]]$x[trainrows,], 
                                               sim$datasets[[1]]$y[trainrows], 
                                               sim$datasets[[1]]$censor[trainrows]))
      
      #the validation dataset
      local.dataframe.test = data.frame(cbind(sim$datasets[[1]]$x[-trainrows,], 
                                              sim$datasets[[1]]$y[-trainrows],
                                              sim$datasets[[1]]$censor[-trainrows]))
      
      colnames(local.dataframe.train) = c(paste0("x.", 1:p), "y", "censor")
      colnames(local.dataframe.test) = c(paste0("x.", 1:p), "y", "censor")
      
      fit.local = coxph(surv.formula, data = local.dataframe.train, x = TRUE, 
                        control = coxph.control(iter.max = 50))
      Model.local <- list("Cox.local" = fit.local)
      
      PecLocal <- pec(object = Model.local, formula = surv.formula,
                      data = local.dataframe.test, traindata = local.dataframe.train, reference = F)
      PredErrorLocal = c(PredErrorLocal, as.numeric(crps(PecLocal))[1])
      
      sim$datasets[[1]] <- list(x = as.matrix(local.dataframe.train[, 1:p]), 
                                y = local.dataframe.train[, (p + 1)],
                                censor = local.dataframe.train[, (p + 2)])
      
      fit.pars = NULL
      
      for(method in methods){
        fit.method.pec = fit.local #we want to reuse the original coxph object, replacing coefficients
        
        fit = joint_fitting(datasets = sim$datasets, family = "cox", method = method, trace = trace)
        fit.pars = rbind(fit.pars, fit$par)
        
        fit.method.pec$coefficients = fit$par #replacing the coefficients
        
        #creating a new pec object
        Modelsfit <- list("Cox.fit" = fit.method.pec)
        PecWeighted <- pec(object = Modelsfit, formula = surv.formula,
                           data = local.dataframe.test, traindata = local.dataframe.train, 
                           reference = F)
        
        if(method == "opt"){ #different fit objects for 
          fit.method.opt = fit.method.pec
        } else if(method == "testingopt"){
          fit.method.testingopt = fit.method.pec
        } 
        
        PredErrorweighted.method = c(PredErrorweighted.method, as.numeric(crps(PecWeighted))[1])
        EstErrorweighted.method = c(EstErrorweighted.method, sum((fit$par - sim$par) ^ 2))
      }
      
      PredErrorWeighted = rbind(PredErrorWeighted, PredErrorweighted.method) #attach est and pred error
      EstErrorLocal = c(EstErrorLocal, sum((fit$par.localonly - sim$par) ^ 2))
      EstErrorWeighted = rbind(EstErrorWeighted, EstErrorweighted.method)
    }
    
    PredErrorRatio <- log(PredErrorLocal/PredErrorWeighted) #matrices of log pred and est error
    EstErrorRatio <- log(EstErrorLocal/EstErrorWeighted)
    
    colnames(PredErrorRatio) = paste("w", methods, sep=".")
    colnames(EstErrorRatio) = paste("w", methods, sep=".") 
    
    return(list(PredError = zapsmall(PredErrorRatio), #return list of error improvements by method
                EstError = zapsmall(EstErrorRatio)))
  }
}