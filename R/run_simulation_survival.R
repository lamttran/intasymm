#' Simulations to assess the relative performance of the integration methods in reducing
#' estimation and prediction error
#'
#' @import pec
#' @import glmnet
#' @import survival
#' @import prodlim
#' @param K The total number of datasets (including the local dataset, which is assumed to be the first dataset)
#' @param p The number of parameters; for now hardcoded to 5 due to requirements of pec() function
#' @param n The sample size of the datasets, should be of length K
#' @param K0 The number of datasets to integrate wrt the local dataset (aka have beta = beta0)
#' @param sigma.beta Additional noise to add to non-integrated or fractionally integrated datasets
#' @param mu.x Mean of the design matrix
#' @param sigma.x Standard deviation of the design matrix
#' @param numfrac Of the K0 to-be-integrated datasets, how many should be fractionally integrated 
#' (i.e. beta has some noise)
#' @param sigma.y Noise to add to linear predictor in calculating exponentially distributed survival times
#' @param seeds Seeds to use to replicate simulations
#' @param methods Methods to run simulations for; refer to joint_fitting_survival documentation for further
#' information
#' @param correlated Designates whether the design matrix columns should be correlated (common in genomic data).
#' @param equalcor Only applies if correlated = TRUE. Designates if the correlations the same in all datasets.
#' @param trainingsize The size of the local dataset (sampled from the first dataset) to be
#' used as training data 
#' @return A matrix of improvements in estimation and prediction error, each column representing a specified
#' method. Values are presented as the log of the local-only error (no integration) divided by the post-
#' integration error; that is positive values represent improvements by performing integration. 
#' @export

run_simulation_survival = function (K = 3, p = 5, n = c(100, rep(200, K - 1)), 
                           K0 = 1, sigma.beta = 1, mu.x = 0, sigma.x = 1, sigma.y = 0, 
                           seeds = c(1:6), numfrac = 0, trainingsize = 50,
                           methods = c("opt", "testingopt", "lasso", "ridge"), 
                           correlated = F, equalcor = F, trace = F) { 
  
  PredErrorLocal = NULL #intialize estimation and prediction error matrices
  PredErrorWeighted = NULL
  EstErrorLocal = NULL
  EstErrorWeighted = NULL
  
  for(i in 1:length(seeds)){
    print(i)
    
    set.seed(seeds[i])
    PredErrorweighted.method = NULL #that iteration's prediction and estimation error matrices
    EstErrorweighted.method = NULL
    
    sim = simulate_survival(K = K, p = p, K0 = K0, n = n, numfrac = numfrac, 
                            sigma.beta = sigma.beta, sigma.y = sigma.y, correlated = correlated, 
                            equalcor = equalcor) #the simulated dataset
    
    #select rows for training data i.e. the "local" dataset
    trainrows = c(sort(sample(1:(nrow(sim$datasets[[1]]$x) - 1), size = trainingsize, 
                              replace = FALSE)), nrow(sim$datasets[[1]]$x)) 
    
    #the local dataset
    local.dataframe.train = data.frame(cbind(sim$datasets[[1]]$x[trainrows,], sim$datasets[[1]]$y[trainrows], 
                                             sim$datasets[[1]]$censor[trainrows]))
    
    #the validation dataset
    local.dataframe.test = data.frame(cbind(sim$datasets[[1]]$x[-trainrows,], sim$datasets[[1]]$y[-trainrows],
                                            sim$datasets[[1]]$censor[-trainrows]))
    
    #setting column names. hard-coded as 5 parameters for now due to pec issues
    colnames(local.dataframe.train) = c("x.1", "x.2", "x.3", "x.4", "x.5", "y", "censor")
    colnames(local.dataframe.test) = c("x.1", "x.2", "x.3", "x.4", "x.5", "y", "censor")
    
    #the local-only model, in this format for pec
    fit.local = coxph(Surv(y, censor) ~ x.1 + x.2 + x.3 + x.4 + x.5,
                      data = local.dataframe.train, x = TRUE, control = coxph.control(iter.max = 50))
    Model.local <- list("Cox.local" = fit.local)
    
    #the local-only pec object and prediction error
    PecLocal <- pec(object = Model.local,
                    formula = Surv(y, censor) ~ x.1 + x.2 + x.3 + x.4 + x.5,
                    data = local.dataframe.test, traindata = local.dataframe.train, reference = F)
    PredErrorLocal = c(PredErrorLocal, as.numeric(crps(PecLocal))[1])
    
    #sim$1datasets is now just the training data
    sim$datasets[[1]] <- list(x = as.matrix(local.dataframe.train[, 1:5]), y = local.dataframe.train[, 6],
                              censor = local.dataframe.train[, 7])
    fit.pars = NULL
    
    for(method in methods){
      fit.method.pec = fit.local #we want to reuse the original coxph object, replacing coefficients
      
      fit = joint_fitting_survival(datasets = sim$datasets, fit_weighted_survival, loocv_weighted_survival, 
                          method = method, trace = T) 
      fit.pars = rbind(fit.pars, fit$par)
      
      fit.method.pec$coefficients = fit$par #replacing the coefficients
      
      #creating a new pec object
      Modelsfit <- list("Cox.fit" = fit.method.pec)
      PecWeighted <- pec(object = Modelsfit,
                         formula = Surv(y, censor) ~ x.1 + x.2 + x.3 + x.4 + x.5,
                         data = local.dataframe.test, traindata = local.dataframe.train, 
                         reference = F)
      
      if(method == "opt"){ #different fit objects for 
        fit.method.opt = fit.method.pec
      } else if(method == "testingopt"){
        fit.method.testingopt = fit.method.pec
      } else if (method == "lasso"){
        fit.method.lasso = fit.method.pec
      } else if (method == "ridge"){
        fit.method.ridge = fit.method.pec
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