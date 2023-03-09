#Code for the imputation of HeiNet with ERGM and for analysis with SAOM and with STERGM

#An equivalent was applied for the stationary SAOM,  Bayesian ERGM and SAOM internal imputation.
#The last one means that the first wave is imputed by null imputation and the later waves by the SAOM internal mechanism.
#In any case the waves two and three are imputed by the SAOM internal mechanism.

#Multiple imputation based on an ERGM to impute the first wave was first applied by Hipp et al. 2015:
#https://www.sciencedirect.com/science/article/abs/pii/S0378873315000027
#This script loosely follows their description

#Please note: The ERGM version 4.1.2 (newest when this analysis was conducted) had a bug that does not allow to reasonably draw simulatiomns with 
#dayadwise restrictions (keep observed) that are necessary for imputation, so Version 3.8.0 was used, but 3.8.0 is NOT compatible with Bergm,
#Therefore 4.1.2 is used for Bergm, please re-install ERGM package if formerly downgraded to 3.8.0 before running Bergm.

#get packages
library(tidyverse)
library(knitr)
library(RSiena)
library(stargazer)
library(ergm) #version 3.8.0
library(tergm)
library(latticeExtra)



#get the data
load("mydata.RData")
#Bekannte_full_w1_ma, Bekannte_full_w2_ma and Bekannte_full_w3_ma are adjecency matrices in the matrix format.
#They have a size of 42x42 and represent the the waves of relational data collected, missings are NA.
#The data also includes attribute data about the tutorial group in statistics each individual participated in in the data.table Attributes_K1.
#Missings in the attributes are not imputed.

# Prepare the imputation of the first wave
wave1imp_ergm <- list() #create list to store the imputed networks
D <- 50 #number of imputations
ergmImputation_numImpEdges <- list() #list for numbers of edges of the imputed networks

#estimate the imputation model
ergmImputation_w1g <- ergm(Bekannte_full_net_w1 ~ edges
                           + mutual
                           + gwesp(decay = log(2), fixed = TRUE)
                           + edgecov(Bekannte_full_w2_maI)
                           + odegree1.5()
                           + nodematch('Statistik'),
                           control = control.ergm(    seed        = 123,
                                                      MCMLE.maxit = 60,
                                                      parallel    = 5,
                                                      CD.maxit    = 15,
                                                      MCMC.samplesize = 4000*5, #multiple of number of chains used
                                                      obs.MCMC.samplesize = 4000*5*4, #should be 4*MCMC.samplesize to make concvergence check work adequately
                                                      MCMC.burnin = 30000,
                                                      MCMC.interval = 4000*5), 
                           verbose = FALSE)

#check AIC and BIC to find a good model
summary(ergmImputation_w1g)

#MCMC diagnostics
mcmc.diagnostics(ergmImputation_w1g)

#Goodness of Fit
ergmImputation_w1g_gof <- gof(ergmImputation_w1g, GOF=~ idegree + odegree + espartners + distance + model)
plot(ergmImputation_w1g_gof)
ergmImputation_w1g_gof
#make sure the p-values are >0


##### Simulation #####

#Hipp et al. (2015): all simulations are based on the same model estimation as it would be too computationally extensive otherwise

ergm_imputation <- simulate(ergmImputation_w1g, nsim=D, constraints = ~observed)
summary(ergm_imputation[[1]])

##### examine the model #####
#t value is supposed to be lower than 0.1 (like SAOM) and tconv.max should be lower than 0.25 (ideally) and definitely lower than 0.35


tk.list <- list()
tconv.max.list  <- list()

for (i in 1:5){ #5 times becaused 5 chains were used
  a <- as.matrix(ergmImputation_w1g$sample[[i]])
  b <- as.matrix(ergmImputation_w1g$sample.obs[[i]])
  sf <- a-b
  
  tk <- abs(apply(sf,2,mean))/apply(a,2,sd)
  tk.list[[i]] <- tk
  tconv.max <- sqrt(t(apply(sf,2,mean) %*% solve(as.matrix(cov(sf))) %*% apply(sf,2,mean)))
  tconv.max.list[[i]] <- tconv.max
  
}

#View results
tk.list
tconv.max.list

tconv.max <- sqrt(t(apply(as.matrix(sf_list_ergm[[1]]),2,mean) %*% solve(as.matrix(cov(as.matrix(sf_list_ergm[[1]])))) %*% apply(as.matrix(sf_list_ergm[[1]]),2,mean)))
tconv.max

##### save the imputed networks #####
for (i in 1:D) {
  wave1imp_ergm[[i]] <- as.matrix.network(ergm_imputation[[i]])
  ergmImputation_numImpEdges[i] <- table(wave1imp_ergm[[i]])[2] #save only number of 1s
}

View(wave1imp_ergm[[1]])
ergmImputation_numImpEdges #number of existing edges in each imputation

save(ergmImputation_w1g, ergmImputation_w1g_gof, ergm_imputation,
     ergmImputation_numImpEdges, wave1imp_ergm, file = "ergm_Imputation_final.RData")



##### Later waves imputed by SAOM internal mechanism #####

#Define Function that will hopefully ensure convergence in RSiena
siena07ToConvergence <- function(alg, dat, eff, ans0 = NULL, threshold = 0.25,
                                 nodes = 1, ...){
  # parameters are:
  # alg, dat, eff: Arguments for siena07: algorithm, data, effects object.
  # ans0: previous answer, if available; used as prevAns in siena07.
  # threshold: largest satisfactory value
  #            for overall maximum convergence ratio (indicating convergence).
  # nodes: number of processes for parallel processing.
  if (!is.null(ans0)) {
    alg$nsub = 6
  }
  numr <- 0 # number of repetitions
  ans <- siena07(alg, data = dat, effects = eff, prevAns = ans0,
                 nbrNodes = nodes, returnDeps = TRUE,
                 useCluster = (nodes >= 2), ...) # the first run
  repeat {
    save(ans, file = paste("ans",numr,".RData",sep = "")) # to be safe
    numr <- numr + 1         # count number of repeated runs
    tm <- ans$tconv.max      # convergence indicator
    cat(numr, tm,"\n")       # report how far we are
    if (tm < threshold) {break}   # success
    if (tm > 10) {break}     # divergence without much hope
    # of returning to good parameter values
    if (numr > 20) {break}  # now it has lasted too long
    alg$nsub <- 1
    alg$n2start <- 2 * (sum(eff$include) + 7) * 2.52**4
    alg$n3 <- 3000
    if (numr == 1) {alg$firstg <- alg$firstg/5}
    ans <- siena07(alg, data = dat, effects = eff, prevAns = ans,
                   nbrNodes = nodes, returnDeps = TRUE,
                   useCluster = (nodes >= 2),...)
  }
  if (tm > threshold)
  {
    stop("Convergence inadequate.\n")
  }
  ans
}

#Since version 1.2-12, Maximum Likelihood (ML) estimation by Rsiena with returnDeps = TRUE returns an edgelist of the final network at the end of the phase 3 simulation.
#The following function, getNet(), uses this edgelist to impute the data.

getNet <- function(observedNet,edgeList) {
  # observedNet = observed network as adjacency with missing data
  # edgeList = edgelist that is returned by siena07(...)$sims
  observedNet[is.na(observedNet)] <- 0
  for (i in 1:nrow(edgeList)) {
    observedNet[edgeList[i,1],edgeList[i,2]] <- 1
  }
  return(observedNet)
}

D <- 50 # number of imputations

#create lists for the imputed networks of wave 2 and 3
wave2imp_ergm <- list()
wave3imp_ergm <- list()

#prepare for SAOM
Statistik_v <- as.vector(Attributes_K1$Statistik)
Statistik  <- coCovar(as.vector(Statistik_v),centered=FALSE)

#set estimation options (same as used for Stationary SAOM and the 2 latter waves after Stat. SAOM)
estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE,
                                           seed = 2214,
                                           n3 = 1000, maxlike = FALSE,
                                           cond = FALSE) #n3 =number of datasets simulated

#set imputation options (same as used for Stationary SAOM and the 2 latter waves after Stat. SAOM)
imputation.options <- sienaAlgorithmCreate(seed = 13848,
                                           useStdInits = FALSE, 
                                           maxlike = TRUE,
                                           cond = FALSE, 
                                           nsub = 0,
                                           simOnly = TRUE,
                                           n3 = 10)

#For each wave estimation is done WIth Method of Moments (MoM) and imputation with Maximum Likelihood (ML)
#Imputation is repeated as often as the first wave had been imputed (D times)

set.seed(1307)
for (i in 1:D) {
  
  cat('imputation',i,'\n')
  
  # now impute wave2
  
  Bekannte_full <- sienaDependent(array(c(wave1imp_ergm[[i]],Bekannte_full_w2_ma),
                                        dim = c(42,42,2)))
  
  Data.w2  <- sienaDataCreate(Bekannte_full, Statistik)
  
  
  effects.twoWaves <- getEffects(Data.w2)
  effects.twoWaves <- includeEffects(effects.twoWaves, 
                                     gwespFF, recip, density, outAct)
  effects.twoWaves <- includeEffects(effects.twoWaves, sameX,
                                     interaction1 = "Statistik")
  
  if (i == 1) {
    period1saom <- siena07ToConvergence(alg = estimation.options,
                                        dat = Data.w2,
                                        eff = effects.twoWaves,
                                        threshold = 0.25)
  } else {
    period1saom <- siena07ToConvergence(alg = estimation.options,
                                        dat = Data.w2,
                                        eff = effects.twoWaves,
                                        threshold = 0.25,
                                        ans0 = period1saom)
  }
  
  sims <- siena07(imputation.options, data = Data.w2,
                  effects = effects.twoWaves,
                  prevAns = period1saom,
                  returnDeps = TRUE)$sims[[10]][[1]]
  
  wave2imp_ergm[[i]] <- getNet(Bekannte_full_w2_ma,sims)
  
  # impute wave 3
  
  Bekannte_full <- sienaDependent(array( c(wave2imp_ergm[[i]], Bekannte_full_w3_ma),
                                         dim = c( 42, 42, 2)))

  Data.w3  <- sienaDataCreate(Bekannte_full, Statistik)
  
  if (i == 1) {
    period2saom <- siena07ToConvergence(alg = estimation.options,
                                        dat = Data.w3,
                                        eff = effects.twoWaves,
                                        threshold = 0.25)
  } else {
    period2saom <- siena07ToConvergence(alg = estimation.options,
                                        dat = Data.w3,
                                        eff = effects.twoWaves,
                                        threshold = 0.25,
                                        ans0 = period2saom)
  }
  
  
  sims <- siena07(imputation.options, data = Data.w3,
                  effects = effects.twoWaves,
                  prevAns = period2saom,
                  returnDeps = TRUE)$sims[[10]][[1]]
  
  wave3imp_ergm[[i]] = getNet(Bekannte_full_w3_ma,sims)
}

#check results and rename to identify results later:

period1saom_ergm <- period1saom
period2saom_ergm <- period2saom

#ckeck covariance
summary(period1saom_ergm)
summary(period2saom_ergm)

#check for autocorrelation
period1saom_ergm$ac
period2saom_ergm$ac

#Goodness of Fit
ergm.saom.model.results.gof.p1 <- list()
ergm.saom.model.results.gof.p1[[1]] <- plot(sienaGOF(period1saom_ergm, varName="Bekannte_full", OutdegreeDistribution))
ergm.saom.model.results.gof.p1[[2]] <- plot(sienaGOF(period1saom_ergm, varName="Bekannte_full", IndegreeDistribution))
#View results
ergm.saom.model.results.gof.p1[[1]]
ergm.saom.model.results.gof.p1[[2]]

ergm.saom.model.results.gof.p2 <- list()
ergm.saom.model.results.gof.p2[[1]] <- plot(sienaGOF(period2saom_ergm, varName="Bekannte_full", OutdegreeDistribution))
ergm.saom.model.results.gof.p2[[2]] <- plot(sienaGOF(period2saom_ergm, varName="Bekannte_full", IndegreeDistribution))
#View results
ergm.saom.model.results.gof.p2[[1]]
ergm.saom.model.results.gof.p2[[2]]


#call a certain imputation, e.g., the 4th imputation of wave 2:
View(wave2imp_ergm[[4]])

#save imputations
save(ergmImputation_w1g, ergmImputation_w1g_gof, ergm_imputation,
     ergmImputation_numImpEdges, wave1imp_ergm, wave2imp_ergm, wave3imp_ergm,
     ergm.saom.model.results.gof.p1, ergm.saom.model.results.gof.p2,
     period1saom_ergm, period2saom_ergm, significance_p1_ergm, significance_p2_ergm,
     file = "ergm_Imputation_final_all_waves.RData")

########################################################################################################################
####################################### Estimating the analysis models #################################################
########################################################################################################################

############### SAOM Analysis for the Stationary ERGM imputation ##################

#run our model on the D complete data sets in a loop
#output is saved in a list
saomResults_ergm <- list()

for (i in 1:D) {
  cat('Imputation',i,'\n')
  
  Bekannte_full <- sienaDependent(array(c(wave1imp_ergm[[i]], wave2imp_ergm[[i]],
                                          wave3imp_ergm[[i]]), dim = c(42, 42, 3)))
  
  Data  <- sienaDataCreate(Bekannte_full, Statistik)
  effectsData <- getEffects(Data)
  effectsData <- includeEffects(effectsData,
                                gwespFF, recip, density, outAct)
  effectsData <- includeEffects(effectsData,  sameX,
                                interaction1 =  "Statistik")
  
  estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 312,
                                             n3 = 1000, maxlike = FALSE, cond = FALSE,
                                             lessMem = TRUE)
  if (i == 1) {
    saomResults_ergm[[i]] <- siena07ToConvergence(alg = estimation.options,
                                                  dat = Data, eff = effectsData,
                                                  threshold = 0.25)
  } else {
    saomResults_ergm[[i]] <- siena07ToConvergence(alg = estimation.options,
                                                  dat = Data, eff = effectsData,
                                                  threshold = 0.25,
                                                  ans0 = saomResults_ergm[[i - 1]])
  }
}

#check rsults of first SAOM
summary(saomResults_ergm[[1]])

##### Combining the Results by Rubin's rules #####

#The combined estimate for the parameters is the average of the estimates of the D analyses

rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}


npar <- sum(effectsData$include) #amount of parameters included

#create an dataframe with all estimated parameters and standard errors
MIResults <- as.data.frame(matrix(,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- saomResults_ergm[[i]]$theta
  MIResults[,i * 2] <-  sqrt(diag(saomResults_ergm[[i]]$covtheta))
}

#Now we get the average covariance structure between the parameters

WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + saomResults_ergm[[i]]$covtheta
}

WDMIs <- (1/D) * WDMIs

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure

finalResults_ergm <- as.data.frame(matrix(,npar,2))
names(finalResults_ergm) <- c("combinedEstimate", "combinedSE")
rownames(finalResults_ergm) <- effectsData$effectName[effectsData$include]
finalResults_ergm$combinedEstimate <- rowMeans(MIResults[,seq(1,2*D,2)])
finalResults_ergm$combinedSE <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                                       rowVar(MIResults[,seq(1,2*D,2)]))

#Odds Ratios
finalResults_ergm$OddsRatio <- exp(finalResults_ergm$combinedEstimate)

#rounding 
finalResults_ergm_rounded <- as.data.frame(round(finalResults_ergm, 3))
View(finalResults_ergm_rounded)

#save results
save(ergmImputation_w1g, ergmImputation_w1g_gof, ergm_imputation,
     ergmImputation_numImpEdges, wave1imp_ergm, wave2imp_ergm, wave3imp_ergm,
     ergm.saom.model.results.gof.p1, ergm.saom.model.results.gof.p2,
     period1saom_ergm, period2saom_ergm, significance_p1_ergm, significance_p2_ergm,
     finalResults_ergm, finalResults_ergm_rounded, MIResults,
     saomResults_ergm, effectsData,
     file = "ergm_Imputation_Analysis_Results.RData")

#get a LaTeX ready table
stargazer(finalResults_ergm_rounded, summary = FALSE, title="SAOM Analyse, 1. Welle imputiert mit ERGM",
          rownames = TRUE)

#significance levels
finalResults_ergm_rounded$sig <- NA
finalResults_ergm_rounded$sig[(abs(finalResults_ergm_rounded$combinedEstimate) - abs(finalResults_ergm_rounded$combinedSE)*1.96)>0]<- "*"
finalResults_ergm_rounded$sig[(abs(finalResults_ergm_rounded$combinedEstimate) - abs(finalResults_ergm_rounded$combinedSE)*2.58)>0]<- "**"
finalResults_ergm_rounded$sig[(abs(finalResults_ergm_rounded$combinedEstimate) - abs(finalResults_ergm_rounded$combinedSE)*3.29)>0]<- "***"

############## Analysis: TERGM Estimation for Stationary SAOM imputation #############

#The analysis is rerun using a separable temporal exponential random graph model. STERGM models formation and dissolution effects separately, hence "separable TERGM".
#here dissolution effects are only included to improve the fit.
#As a model for all data at once did not converge, the model is estimated for both periods separately.

#load the imputed data sets if they are not already in the environment
load("ergm_Imputation_final_all_waves.RData")

#number of imputations (skip if initialized above)
D <- 50

#Convert the SAOM simulations to network for tergm to use
wave1imp_net_ergm <- list()
wave2imp_net_ergm <- list()
wave3imp_net_ergm <- list()

for (i in 1:D) {
  wave1imp_net_ergm[[i]] <- as.network(wave1imp_ergm[[i]])
  wave2imp_net_ergm[[i]] <- as.network(wave2imp_ergm[[i]])
  wave3imp_net_ergm[[i]] <- as.network(wave3imp_ergm[[i]])
}

#Add statistics tutorial attribute to all the simulated networks
set.vertex.attribute(Bekannte_full_net_w1, 'Statistik', Attributes_K1$Statistik_Zero)

for (i in 1:D) {
  set.vertex.attribute(wave1imp_net_ergm[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
  set.vertex.attribute(wave2imp_net_ergm[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
  set.vertex.attribute(wave3imp_net_ergm[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
}


#period 1: wave 1 and 2
Result_tergm_ergm_w12 <- list()

set.seed(456)
for (i in 1:D) {
  cat('Imputation',i,'\n')
  Result_tergm_ergm_w12[[i]] <- tergm(list(wave1imp_net_ergm[[i]], wave2imp_net_ergm[[i]])~
                                        Form(~ edges
                                             + mutual
                                             + gwesp(decay = log(2), fixed = TRUE)
                                             + odegree1.5()
                                             + nodematch('Statistik'))+
                                        Diss(~ edges
                                             + mutual
                                             + nodematch('Statistik')),
                                      estimate="CMLE")
}
#see result with 1st imputation
summary(Result_tergm_ergm_w12[[1]])


#period 2: wave 2 and 3
Result_tergm_ergm_w23 <- list()

set.seed(456)
for (i in 1:D) {
  cat('Imputation',i,'\n')
  Result_tergm_ergm_w23[[i]] <- tergm(list(wave2imp_net_ergm[[i]], wave3imp_net_ergm[[i]])~
                                        Form(~ edges
                                             + mutual
                                             + gwesp(decay = log(2), fixed = TRUE)
                                             + odegree1.5()
                                             + nodematch('Statistik'))+
                                        Diss(~ edges
                                             + mutual
                                             + nodematch('Statistik')),
                                      estimate="CMLE")
}
#see result with 1st imputation
summary(Result_tergm_ergm_w23[[1]])

##### create data set for pooling wave 1 and 2 #####

npar <- 8 #number of parameters included in the TERGM
TERGM_Results_ergm_combined_w12 <- as.data.frame(matrix(,npar,(4 * D))) #data frame for results

for (i in 1:D) {
  names(TERGM_Results_ergm_combined_w12)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(TERGM_Results_ergm_combined_w12)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  TERGM_Results_ergm_combined_w12[,i * 2 - 1] <- summary(Result_tergm_ergm_w12[[i]])$coefs$"Estimate"
  TERGM_Results_ergm_combined_w12[,i * 2] <-  summary(Result_tergm_ergm_w12[[i]])$coefs$"Std. Error"
}

#Now we get the average covariance structure between the parameters
WDMIs_w12 <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs_w12 <- WDMIs_w12 + Result_tergm_ergm_w12[[i]]$est.cov
}

WDMIs_w12 <- (1/D) * WDMIs_w12

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

finalResults_TERGM_ergm_w12 <- as.data.frame(matrix(,npar,2))
names(finalResults_TERGM_ergm_w12) <- c("combinedEstimate", "combinedSE")
rownames(finalResults_TERGM_ergm_w12) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                           "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_ergm_w12$combinedEstimate <- rowMeans(TERGM_Results_ergm_combined_w12[,seq(1,2*D,2)])
finalResults_TERGM_ergm_w12$combinedSE <- sqrt(diag(WDMIs_w12) + ((D + 1)/D) *
                                                 rowVar(TERGM_Results_ergm_combined_w12[,seq(1,2*D,2)]))

#Odds Ratios
finalResults_TERGM_ergm_w12$OddsRatio <- exp(finalResults_TERGM_ergm_w12$combinedEstimate)

#rounding
finalResults_TERGM_ergm_rounded_w12 <- as.data.frame(round(finalResults_TERGM_ergm_w12, 3))
View(finalResults_TERGM_ergm_rounded_w12)

##### create data set for pooling wave 2 and 3 #####
npar <- 8 #number of parameters included in the TERGM
TERGM_Results_ergm_combined_w23 <- as.data.frame(matrix(,npar,(4 * D))) #data frame for results

for (i in 1:D) {
  names(TERGM_Results_ergm_combined_w23)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(TERGM_Results_ergm_combined_w23)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  TERGM_Results_ergm_combined_w23[,i * 2 - 1] <- summary(Result_tergm_ergm_w23[[i]])$coefs$"Estimate"
  TERGM_Results_ergm_combined_w23[,i * 2] <-  summary(Result_tergm_ergm_w23[[i]])$coefs$"Std. Error"
}

#Now we get the average covariance structure between the parameters

WDMIs_w23 <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs_w23 <- WDMIs_w23 + Result_tergm_ergm_w23[[i]]$est.cov
}

WDMIs_w23 <- (1/D) * WDMIs_w23

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure
finalResults_TERGM_ergm_w23 <- as.data.frame(matrix(,npar,2))
names(finalResults_TERGM_ergm_w23) <- c("combinedEstimate", "combinedSE")
rownames(finalResults_TERGM_ergm_w23) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                           "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_ergm_w23$combinedEstimate <- rowMeans(TERGM_Results_ergm_combined_w23[,seq(1,2*D,2)])
finalResults_TERGM_ergm_w23$combinedSE <- sqrt(diag(WDMIs_w23) + ((D + 1)/D) *
                                                 rowVar(TERGM_Results_ergm_combined_w23[,seq(1,2*D,2)]))

#Odds Ratios
finalResults_TERGM_ergm_w23$OddsRatio <- exp(finalResults_TERGM_ergm_w23$combinedEstimate)

#rounding
finalResults_TERGM_ergm_rounded_w23 <- as.data.frame(round(finalResults_TERGM_ergm_w23, 3))
View(finalResults_TERGM_ergm_rounded_w23)

save(wave1imp_net_ergm, wave2imp_net_ergm, wave3imp_net_ergm,
     finalResults_TERGM_ergm_w12, finalResults_TERGM_ergm_rounded_w12, TERGM_Results_ergm_combined_w12,
     Result_tergm_ergm_w12, finalResults_TERGM_ergm_w23, finalResults_TERGM_ergm_rounded_w23,
     TERGM_Results_ergm_combined_w23, Result_tergm_ergm_w23,
     file = "ERGM_TERGM_Analysis_Results.RData")

###### Assessing goodness of fit #####
gof_tergm_ergm_w12 <- gof(Result_tergm_ergm_w12[[6]])
gof_tergm_ergm_w12 
plot(gof_tergm_ergm_w12)

mcmc.diagnostics(Result_tergm_ergm_w12[[6]])
mcmc.diagnostics(Result_tergm_ergm_w12[[4]])
mcmc.diagnostics(Result_tergm_ergm_w23[[6]])
mcmc.diagnostics(Result_tergm_ergm_w23[[4]])

gof_tergm_ergm_w23 <- gof(Result_tergm_ergm_w23[[6]])
gof_tergm_ergm_w23 
plot(gof_tergm_ergm_w23)

#make LaTeX ready table
rownames(finalResults_TERGM_ergm_rounded_w12) <- c("Form edges", "Form mutual", "Form gwesp.fixed.0.69", "Form odegree1.5", "Form nodematch.Statistik",
                                                   "Diss edges", "Diss mutual", "Diss nodematch.Statistik")

#combined
finalResults_TERGM_ergm_rounded_w12_and_w23 <- finalResults_TERGM_ergm_rounded_w12
finalResults_TERGM_ergm_rounded_w12_and_w23$OddsRatio12 <- round(exp(finalResults_TERGM_ergm_rounded_w12$combinedEstimate), 3)
finalResults_TERGM_ergm_rounded_w12_and_w23$combinedEstimate23 <- finalResults_TERGM_ergm_rounded_w23$combinedEstimate
finalResults_TERGM_ergm_rounded_w12_and_w23$combinedSE23 <- finalResults_TERGM_ergm_rounded_w23$combinedSE
finalResults_TERGM_ergm_rounded_w12_and_w23$OddsRatio23 <- round(exp(finalResults_TERGM_ergm_rounded_w23$combinedEstimate), 3)

stargazer(finalResults_TERGM_ergm_rounded_w12_and_w23, summary = FALSE, title="TERGM Anayse, 1st wave imputed by ERGM",
          rownames = TRUE, column.labels= c("Period 1", "Period 2"), column.separate = c(3,3))


#significance levels

finalResults_TERGM_ergm_rounded_w12_and_w23$sig_w12 <- NA
finalResults_TERGM_ergm_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedSE)*1.96)>0]<- "*"
finalResults_TERGM_ergm_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedSE)*2.58)>0]<- "**"
finalResults_TERGM_ergm_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedSE)*3.29)>0]<- "***"

finalResults_TERGM_ergm_rounded_w12_and_w23$sig_w23 <- NA
finalResults_TERGM_ergm_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedSE23)*1.96)>0]<- "*"
finalResults_TERGM_ergm_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedSE23)*2.58)>0]<- "**"
finalResults_TERGM_ergm_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_ergm_rounded_w12_and_w23$combinedSE23)*3.29)>0]<- "***"



