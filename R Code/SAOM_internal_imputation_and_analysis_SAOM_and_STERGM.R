#Code for the imputation of HeiNet with the SAOM internal mechanism and for analysis with SAOM and with STERGM

#An equivalent was applied for the stationary SAOM, ERGM and SAOM internal imputation.
#The last one means that the first wave is imputed by null imputation and the later waves by the SAOM internal mechanism.
#In any case the waves two and three are imputed by the SAOM internal mechanism. In this script the fist wave is null imputed (NA == 0)

#get packages
library(tidyverse)
library(knitr)
library(RSiena)
library(stargazer)
library(ergm)
library(tergm)

#get the data
load("mydata.RData")
#Bekannte_full_w1_ma, Bekannte_full_w2_ma and Bekannte_full_w3_ma are adjecency matrices in the matrix format.
#They have a size of 42x42 and represent the the waves of relational data collected, missings are NA.
#The data also includes attribute data about the tutorial group in statistics each individual participated in in the data.table Attributes_K1.
#Missings in the attributes are not imputed.

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

D <- 50 #number of imputations

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

#create lists for the imputed networks of wave 1, 2 and 3
wave2imp_SAOM_Int <- list()
wave3imp_SAOM_Int <- list()

#prepare for SAOM
Statistik_v <- as.vector(Attributes_K1$Statistik)
Statistik  <- coCovar(as.vector(Statistik_v),centered=FALSE)

#For each wave estimation is done WIth Method of Moments (MoM) and imputation with Maximum Likelihood (ML)

set.seed(1307)
for (i in 1:D) {
  
  cat('imputation',i,'\n')
  
  # now impute wave2
  
  Bekannte_full <- sienaDependent(array(c(Bekannte_full_w1_ma, Bekannte_full_w2_ma), #wave 1 just kept with missings
                                        dim = c(42,42,2)))
  Data.w2  <- sienaDataCreate(Bekannte_full, Statistik)
  
  
  effects.twoWaves <- getEffects(Data.w2)
  effects.twoWaves <- includeEffects(effects.twoWaves,
                                     gwespFF, recip, density, outAct)
  effects.twoWaves <- includeEffects(effects.twoWaves,  sameX,
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
  
  wave2imp_SAOM_Int[[i]] <- getNet(Bekannte_full_w2_ma,sims)
  
  # impute wave 3
  
  Bekannte_full <- sienaDependent(array( c(wave2imp_SAOM_Int[[i]], Bekannte_full_w3_ma),
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
  
  wave3imp_SAOM_Int[[i]] = getNet(Bekannte_full_w3_ma,sims)
}

#check results:

period1saom_SAOM_Int <- period1saom
period2saom_SAOM_Int <- period2saom

#check for collinearity
summary(period1saom_SAOM_Int) 
summary(period2saom_SAOM_Int) 

#check for autocorrelation
period1saom_SAOM_Int$ac
period2saom_SAOM_Int$ac

#Goodness of Fit
saom_SAOM_Int.model.results.gof.p1 <- list()
saom_SAOM_Int.model.results.gof.p1[[1]] <- plot(sienaGOF(period1saom_SAOM_Int, varName="Bekannte_full", OutdegreeDistribution))
saom_SAOM_Int.model.results.gof.p1[[2]] <- plot(sienaGOF(period1saom_SAOM_Int, varName="Bekannte_full", IndegreeDistribution))
#View results
saom_SAOM_Int.model.results.gof.p1[[1]]
saom_SAOM_Int.model.results.gof.p1[[2]]

saom_SAOM_Int.model.results.gof.p2 <- list()
saom_SAOM_Int.model.results.gof.p2[[1]] <- plot(sienaGOF(period2saom_SAOM_Int, varName="Bekannte_full", OutdegreeDistribution))
saom_SAOM_Int.model.results.gof.p2[[2]] <- plot(sienaGOF(period2saom_SAOM_Int, varName="Bekannte_full", IndegreeDistribution))
#View results
saom_SAOM_Int.model.results.gof.p2[[1]]
saom_SAOM_Int.model.results.gof.p2[[2]]

#save the model results
save(wave2imp_SAOM_Int, wave3imp_SAOM_Int,
     saom_SAOM_Int.model.results.gof.p1, saom_SAOM_Int.model.results.gof.p2,
     period1saom_SAOM_Int, period2saom_SAOM_Int, file = "SAOM_internal_Imputation.RData")


########################################################################################################################
####################################### Estimating the analysis models #################################################
########################################################################################################################

############### SAOM Analysis for the SAOM internal imputation ##################

#run model on the D complete data sets in a loop
#output is saved in a list
saomResults <- list()

for (i in 1:D) {
  cat('Imputation',i,'\n')
  
  Bekannte_full <- sienaDependent(array(c(Bekannte_full_w1_ma, wave2imp_SAOM_Int[[i]],
                                          wave3imp_SAOM_Int[[i]]), dim = c( 42, 42, 3)))
  
  Data  <- sienaDataCreate(Bekannte_full, Statistik)
  effectsData <- getEffects(Data)
  effectsData <- includeEffects(effectsData,
                                gwespFF, recip, density, outAct)
  effectsData <- includeEffects(effectsData,  sameX,
                                interaction1 =  "Statistik")
  
  estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 312,
                                             n3 = 1000, maxlike = FALSE, cond = FALSE,
                                             lessMem = FALSE)
  if (i == 1) {
    saomResults[[i]] <- siena07ToConvergence(alg = estimation.options,
                                             dat = Data, eff = effectsData,
                                             threshold = 0.25)
  } else {
    saomResults[[i]] <- siena07ToConvergence(alg = estimation.options,
                                             dat = Data, eff = effectsData,
                                             threshold = 0.25,
                                             ans0 = saomResults[[i - 1]])
  }
}

saomResults[[1]]

#colliearity check
summary(saomResults[[1]])

#check for autocorrelation
saomResults[[1]]$ac

#goodness of fit
saom.results.gof.result_SAOM_Int <- list()
saom.results.gof.result_SAOM_Int[[1]] <- plot(sienaGOF(saomResults[[1]], varName="Bekannte_full", OutdegreeDistribution))
saom.results.gof.result_SAOM_Int[[2]] <- plot(sienaGOF(saomResults[[1]], varName="Bekannte_full", IndegreeDistribution))
saom.results.gof.result_SAOM_Int[[1]] 
saom.results.gof.result_SAOM_Int[[2]] 

#rename for identification
saomResults_SAOM_Int <- saomResults

############# Combining the Results by Rubin's rules ################

#The combined estimate for the parameters is the average of the estimates of the D analyses
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

npar <- sum(effectsData$include) #amount of parameters included

#create an dataframe with all estimated parameters and standard errors
MIResults_SAOM_Int <- as.data.frame(matrix(,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults_SAOM_Int)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults_SAOM_Int)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults_SAOM_Int[,i * 2 - 1] <- saomResults_SAOM_Int[[i]]$theta
  MIResults_SAOM_Int[,i * 2] <-  sqrt(diag(saomResults_SAOM_Int[[i]]$covtheta))
}


#Now we get the average covariance structure between the parameters

WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + saomResults_SAOM_Int[[i]]$covtheta
}

WDMIs <- (1/D) * WDMIs

#Using Rubin’s Rules combine the parameters and standard errors and complete the procedure

finalResults_SAOM_Int <- as.data.frame(matrix(,npar,2))
names(finalResults_SAOM_Int) <- c("combinedEstimate", "combinedSE")
rownames(finalResults_SAOM_Int) <- effectsData$effectName[effectsData$include]
finalResults_SAOM_Int$combinedEstimate <- rowMeans(MIResults_SAOM_Int[,seq(1,2*D,2)])
finalResults_SAOM_Int$combinedSE <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                                           rowVar(MIResults_SAOM_Int[,seq(1,2*D,2)]))

#Odds Ratios
finalResults_SAOM_Int$OddsRatio <- exp(finalResults_SAOM_Int$combinedEstimate) #not meaningfull for rate constants

finalResults_SAOM_Int_rounded <- as.data.frame(round(finalResults_SAOM_Int, 3))
View(finalResults_SAOM_Int_rounded)

#save the model results
save(saomResults, finalResults_SAOM_Int, effectsData, MIResults_SAOM_Int, npar, WDMIs,
     finalResults_SAOM_Int_rounded, significance_SAOM_Int_Results, saom.results.gof.result_SAOM_Int, 
     file = "SAOM_Int_Imputation_Analysis_Results.RData")

#rounding
finalResults_SAOM_Int_rounded <- as.data.frame(round(finalResults_SAOM_Int_rounded, 3))
View(finalResults_SAOM_Int_rounded)

#get LaTeX ready table
library(stargazer)
stargazer(finalResults_SAOM_Int_rounded, summary = FALSE, title="SAOM analysis, 1st wave imputed by SAOM internal mechanism",
          rownames = TRUE)

#Significance levels
finalResults_SAOM_Int_rounded$sig <- NA
finalResults_SAOM_Int_rounded$sig[(abs(finalResults_SAOM_Int_rounded$combinedEstimate) - abs(finalResults_SAOM_Int_rounded$combinedSE)*1.96)>0]<- "*"
finalResults_SAOM_Int_rounded$sig[(abs(finalResults_SAOM_Int_rounded$combinedEstimate) - abs(finalResults_SAOM_Int_rounded$combinedSE)*2.58)>0]<- "**"
finalResults_SAOM_Int_rounded$sig[(abs(finalResults_SAOM_Int_rounded$combinedEstimate) - abs(finalResults_SAOM_Int_rounded$combinedSE)*3.29)>0]<- "***"

##### TERGM analysis for SAOM internally imputed networks #####

#The analysis is rerun using a separable temporal exponential random graph model. STERGM models formation and dissolution effects separately, hence "separable TERGM".
#here dissolution effects are only included to improve the fit.
#As a model for all data at once did not converge, the model is estimated for both periods separately.

#load the imputed data sets if they are not already in the environment
load("SAOM_internal_Imputation.RData")

#number of imputations (skip if initialized above)
D <- 50

#Convert the SAOM simulations to network for tergm to use
wave2imp_net_SAOM_Int <- list()
wave3imp_net_SAOM_Int <- list()

for (i in 1:D) {
  wave2imp_net_SAOM_Int[[i]] <- as.network(wave2imp_SAOM_Int[[i]])
  wave3imp_net_SAOM_Int[[i]] <- as.network(wave3imp_SAOM_Int[[i]])
}

#make an NA = 0 version of welle 1
Bekannte_full_w1_ma_null <- Bekannte_full_w1_ma
Bekannte_full_w1_ma_null[is.na(Bekannte_full_w1_ma_null)] <- 0
Bekannte_full_net_w1_null <- as.network(Bekannte_full_w1_ma_null)

#Add Statistik attribute to all the simulated networks
set.vertex.attribute(Bekannte_full_net_w1_null, 'Statistik', Attributes_K1$Statistik_Zero)

for (i in 1:D) {
  set.vertex.attribute(wave2imp_net_SAOM_Int[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
  set.vertex.attribute(wave3imp_net_SAOM_Int[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
}


#Period 1: wave 1 and 2
Result_tergm_SAOM_Int_w12 <- list()

for (i in 1:D) {
  cat('Imputation',i,'\n')
  Result_tergm_SAOM_Int_w12[[i]] <- tergm(list(Bekannte_full_net_w1_null, wave2imp_net_SAOM_Int[[i]])~
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
summary(Result_tergm_SAOM_Int_w12[[1]])

#Period 2:wave 2 and 3
Result_tergm_SAOM_Int_w23 <- list()

for (i in 1:D) {
  cat('Imputation',i,'\n')
  Result_tergm_SAOM_Int_w23[[i]] <- tergm(list(wave2imp_net_SAOM_Int[[i]], wave3imp_net_SAOM_Int[[i]])~
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
summary(Result_tergm_SAOM_Int_w23[[1]])

Result_tergm_SAOM_Int_w12[[4]]$est.cov #covariance matrix of the estimates

####### Make data set for pooling wave 1 and 2 ######
npar <- 8 #number of parameters included in the TERGM
TERGM_Results_SAOM_Int_combined_w12 <- as.data.frame(matrix(,npar,(4 * D))) #data frame for results

for (i in 1:D) {
  names(TERGM_Results_SAOM_Int_combined_w12)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(TERGM_Results_SAOM_Int_combined_w12)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  TERGM_Results_SAOM_Int_combined_w12[,i * 2 - 1] <- summary(Result_tergm_SAOM_Int_w12[[i]])$coefs$"Estimate"
  TERGM_Results_SAOM_Int_combined_w12[,i * 2] <-  summary(Result_tergm_SAOM_Int_w12[[i]])$coefs$"Std. Error"
}

#Now we get the average covariance structure between the parameters

WDMIs_w12 <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs_w12 <- WDMIs_w12 + Result_tergm_SAOM_Int_w12[[i]]$est.cov
}

WDMIs_w12 <- (1/D) * WDMIs_w12

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}


finalResults_TERGM_SAOM_Int_w12 <- as.data.frame(matrix(,npar,2))
names(finalResults_TERGM_SAOM_Int_w12) <- c("combinedEstimate", "combinedSE")
rownames(finalResults_TERGM_SAOM_Int_w12) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                               "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_SAOM_Int_w12$combinedEstimate <- rowMeans(TERGM_Results_SAOM_Int_combined_w12[,seq(1,2*D,2)])
finalResults_TERGM_SAOM_Int_w12$combinedSE <- sqrt(diag(WDMIs_w12) + ((D + 1)/D) *
                                                     rowVar(TERGM_Results_SAOM_Int_combined_w12[,seq(1,2*D,2)]))

#Odds Ratios
finalResults_TERGM_SAOM_Int_w12$OddsRatio <- exp(finalResults_TERGM_SAOM_Int_w12$combinedEstimate)

#rounding
finalResults_TERGM_SAOM_Int_rounded_w12 <- as.data.frame(round(finalResults_TERGM_SAOM_Int_w12, 3))
View(finalResults_TERGM_SAOM_Int_rounded_w12)

####### Make data set for pooling wave 2 and 3 ######
npar <- 8 #number of parameters included in the TERGM
TERGM_Results_SAOM_Int_combined_w23 <- as.data.frame(matrix(,npar,(4 * D))) #data frame for results

for (i in 1:D) {
  names(TERGM_Results_SAOM_Int_combined_w23)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(TERGM_Results_SAOM_Int_combined_w23)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  TERGM_Results_SAOM_Int_combined_w23[,i * 2 - 1] <- summary(Result_tergm_SAOM_Int_w23[[i]])$coefs$"Estimate"
  TERGM_Results_SAOM_Int_combined_w23[,i * 2] <-  summary(Result_tergm_SAOM_Int_w23[[i]])$coefs$"Std. Error"
}

#Now we get the average covariance structure between the parameters

WDMIs_w23 <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs_w23 <- WDMIs_w23 + Result_tergm_SAOM_Int_w23[[i]]$est.cov
}

WDMIs_w23 <- (1/D) * WDMIs_w23

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure
finalResults_TERGM_SAOM_Int_w23 <- as.data.frame(matrix(,npar,2))
names(finalResults_TERGM_SAOM_Int_w23) <- c("combinedEstimate", "combinedSE")
rownames(finalResults_TERGM_SAOM_Int_w23) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                               "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_SAOM_Int_w23$combinedEstimate <- rowMeans(TERGM_Results_SAOM_Int_combined_w23[,seq(1,2*D,2)])
finalResults_TERGM_SAOM_Int_w23$combinedSE <- sqrt(diag(WDMIs_w23) + ((D + 1)/D) *
                                                     rowVar(TERGM_Results_SAOM_Int_combined_w23[,seq(1,2*D,2)]))

#Odds Ratios
finalResults_TERGM_SAOM_Int_w23$OddsRatio <- exp(finalResults_TERGM_SAOM_Int_w23$combinedEstimate)

#rounding
finalResults_TERGM_SAOM_Int_rounded_w23 <- as.data.frame(round(finalResults_TERGM_SAOM_Int_w23, 3))
View(finalResults_TERGM_SAOM_Int_rounded_w23)

#safe results
save(wave2imp_net_SAOM_Int, wave3imp_net_SAOM_Int, Bekannte_full_w1_ma_null, Bekannte_full_net_w1_null,
     finalResults_TERGM_SAOM_Int_w12, finalResults_TERGM_SAOM_Int_rounded_w12, TERGM_Results_SAOM_Int_combined_w12,
     Result_tergm_SAOM_Int_w12, finalResults_TERGM_SAOM_Int_w23, finalResults_TERGM_SAOM_Int_rounded_w23, TERGM_Results_SAOM_Int_combined_w23,
     Result_tergm_SAOM_Int_w23,
     file = "SAOM_Int_TERGM_Analysis_Results.RData")

load("SAOM_Int_TERGM_Analysis_Results.RData")


##### Assessing goodness of fit #####
gof_tergm_saom_int_w12 <- gof(Result_tergm_SAOM_Int_w12[[6]]) #e.g., for imputation no. 6
gof_tergm_saom_int_w12
plot(gof_tergm_saom_int_w12) 

gof_tergm_saom_int_w23 <- gof(Result_tergm_SAOM_Int_w23[[6]])
gof_tergm_saom_int_w23
plot(gof_tergm_saom_int_w23) 

#LaTeX ready table
rownames(finalResults_TERGM_SAOM_Int_rounded_w12) <- c("F edges", "F mutual", "F gwesp(0.69)", "F odegree1.5", "F Fokus Statistik",
                                                       "D edges", "D mutual", "D Fokus Statistik")

#combined
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23 <- finalResults_TERGM_SAOM_Int_rounded_w12
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$OddsRatio <- round(exp(finalResults_TERGM_SAOM_Int_rounded_w12$combinedEstimate), 3)
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedEstimate23 <- finalResults_TERGM_SAOM_Int_rounded_w23$combinedEstimate
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedSE23 <- finalResults_TERGM_SAOM_Int_rounded_w23$combinedSE
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$OddsRatio23 <- round(exp(finalResults_TERGM_SAOM_Int_rounded_w23$combinedEstimate), 3)

stargazer(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23, summary = FALSE, title="TERGM analysis imputed by SAOM internal mechanism",
          rownames = TRUE, column.labels= c("Period 1", "Period 2"), column.separate = c(3,3))


#significance level
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$sig_w12 <- NA
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedSE)*1.96)>0]<- "*"
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedSE)*2.58)>0]<- "**"
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedSE)*3.29)>0]<- "***"

finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$sig_w23 <- NA
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedSE23)*1.96)>0]<- "*"
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedSE23)*2.58)>0]<- "**"
finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_SAOM_Int_rounded_w12_and_w23$combinedSE23)*3.29)>0]<- "***"





