##This script assumes you have loaded the cancersets.RData file into the global environment
##Loading in necessary packages (installed when intasymm is installed)
library(survival)
library(pec)

##Setting a timer
ptm <- proc.time()

##Initializing matrices for error, weights, run-time
CumPredErrorlocalonly = NULL
CumPredErrorweighted = NULL
CumWeights = NULL
CumTime = NULL
fit.methods.all = NULL

##The weighted integration code
for(i in c(1)){ #Our paper used seeds 1:50 for each local data sample size
  set.seed(i)
  print(i)
  ##To predict survival over the testing set, the training data MUST include the latest patient
  ##The sample is training - 1 observations from the 1st to n-1th row of the blca set + the nth row
  ##The paper used size 30 to 100 in increments of 10, so size here is 29, 39 ... 99
  trainrows = c(sort(sample(1:(nrow(blca) - 1), size = 29, 
                            replace = FALSE)), nrow(blca))
  
  ##the new local dataset is the training blca data
  blcatrain = blca[trainrows,]
  blcatest = blca[-trainrows,]
  
  ##Making the dataset as a list of the individual cancer datasets
  datasets = list(list(x = as.matrix(blcatrain[,3:7]), y = blcatrain[,1], censor = blcatrain[,2]), 
                  list(x = as.matrix(brca[,3:7]), y = brca[,1], censor = brca[,2]), 
                  list(x = as.matrix(gbm[,3:7]), y = gbm[,1], censor = gbm[,2]), 
                  list(x = as.matrix(hnsc[,3:7]), y = hnsc[,1], censor = hnsc[,2]), 
                  list(x = as.matrix(laml[,3:7]), y = laml[,1], censor = laml[,2]), 
                  list(x = as.matrix(luad[,3:7]), y = luad[,1], censor = luad[,2]), 
                  list(x = as.matrix(lusc[,3:7]), y = lusc[,1], censor = lusc[,2]), 
                  list(x = as.matrix(ov[,3:7]), y = ov[,1], censor = ov[,2]), 
                  list(x = as.matrix(paad[,3:7]), y = paad[,1], censor = paad[,2]))
  
  ##The local-only model in a pec object to calculate prediction error
  fit.local = coxph(Surv(y, delta) ~ TMX3 + RIBC1 + LOC115581 + HCCAT5 + LOC644151, data = blcatrain, 
                    control = coxph.control(iter.max = 100), x = TRUE)
  
  Model.local <- list("Cox.local" = fit.local)
  
  PredErrorLocal <- pec(object = Model.local,
                        formula = Surv(y, delta) ~ TMX3 + RIBC1 + LOC115581 + HCCAT5 + LOC644151,
                        data = blcatest, traindata = blcatrain, reference = F)
  
  CumPredErrorlocalonly = c(CumPredErrorlocalonly, as.numeric(crps(PredErrorLocal))[1])
  
  ##Code to evaluate prediction error for the opt and testingopt methods
  CumPredErrorweighted.method = NULL
  CumWeights.method = NULL
  CumTime.method = NULL
  fit.pars = NULL
  
  ##For the paper, use both "opt" and "testingopt"
  ##I.e. for(method in c("opt", "testingopt"))
  for(method in c("testingopt")){ 
    ptmmethodstart = proc.time()
    
    ##The post-integration model for opt/testingopt methods
    fit =  joint_fitting_survival(datasets = datasets, fit_weighted_survival, 
                                  loocv_weighted_survival, method = method, 
                                  trace = T)
    
    CumWeights.method = rbind(CumWeights.method, fit$weights)
    
    ##Replacing the coefficient information of the local model with the weighted models
    ##This is because pec requires the usage of a coxph object
    fit.method = fit.local
    fit.pars = rbind(fit.pars, fit$par)
    fit.method$coefficients = fit$par
    Modelsfit <- list("Cox.fit" = fit.method)
    
    ##Calculating prediction error for the weighted integration methods
    PredErrorweighted.method <- pec(object = Modelsfit,
                                    formula = Surv(y, delta) ~ TMX3 + RIBC1 + LOC115581 + HCCAT5 + LOC644151,
                                    data = blcatest, traindata = blcatrain, 
                                    reference = F)
    
    
    CumPredErrorweighted.method = c(CumPredErrorweighted.method, 
                                    as.numeric(crps(PredErrorweighted.method))[1])

    ptmmethodend = proc.time() - ptmmethodstart
    CumTime.method = c(CumTime.method, as.numeric(ptmmethodend[3]))
  }
  
  ##Adding iteration's results to an outputted table
  CumPredErrorweighted = rbind(CumPredErrorweighted, CumPredErrorweighted.method)
  CumWeights = rbind(CumWeights, CumWeights.method)
  CumTime = rbind(CumTime, CumTime.method)
  trace.table = rbind(fit$par.localonly, fit.pars)
  rownames(trace.table) = c("localonly", paste("weighted", method, sep="."))
  fit.methods.all = rbind(fit.methods.all, trace.table)
  print(trace.table)
}

##Ending timer
ptm1 = proc.time() - ptm

##Outputting results of real data analysis
colnames(CumPredErrorweighted) = paste("w", method, sep=".")
boxplot(cbind(CumPredErrorlocalonly, CumPredErrorweighted), 
        names = c("localonly", paste("w", method, sep=".")), 
        main="CumPredError")
print("ASE localonly vs weighted")
print(summary(CumPredErrorlocalonly/CumPredErrorweighted))
print(ptm1)
