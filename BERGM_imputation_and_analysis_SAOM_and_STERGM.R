#Code for the imputation of HeiNet with Bayesian ERGM and for analysis with SAOM and with STERGM

#An equivalent was applied for thestationary SAOM, ERGM and SAOM internal imputation.
#The last one means that the first wave is imputed by null imputation and the later waves by the SAOM internal mechanism.
#In any case the waves two and three are imputed by the SAOM internal mechanism.

#Multiple Imputation for RSiena with a model based approach for the first wave was first applied by Krause et.al.:
#https://www.stats.ox.ac.uk/~snijders/siena/AdSUMMissingDataMD.html
#The imputation and pooling code was used from Krause (2019), with some adjustments

#get packages
library(tidyverse)
library(knitr)
library(RSiena)
library(stargazer)
library(ergm)
library(tergm)
library(Bergm)

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

############# Prepare the imputation of the first wave
wave1imp_bergm <- list() #create list to store the imputed network
D <- 50 #number of imputations

#No centering of covariates in this example because they are factor variables and not metric, otherwise manually center the covariate

#Priors: Rather uninformed priors in this example, setting prior variances
priorVar <- diag(7,6,6) # 7 is the prior variance specified and six is the number of effects, the diagonal is switched to 7 and teh covariance priors remain 0

bergmImputation <- bergmM(Bekannte_full_net_w1 ~ edges
                               + mutual
                               + gwesp(decay = log(2), fixed = TRUE) #OTP
                               + odegree1.5()
                               + edgecov(Bekannte_full_w2_maI)
                               + nodematch('Statistik'),
                               burn.in = 30000, #ERGM Hipp was run with 30000
                               aux.iters = 1000, #n*n with n being number of nodes would be 1764 (see Caimo/Friel 2011 p. 54), but does increase autocorrelation in this case so not done
                               main.iters = D*150,
                               nchains =5, #increasing chains does not help autocorrelation
                               nImp = D,
                               gamma = .3, 
                               prior.sigma = priorVar,
                               seed = 320)

#results
summary(bergmImputation)
plot(bergmImputation, lag = 300) #autocorrelation declines quite quickly


#list of number of imputed edges
bergmImputation_numImpEdges <- list()

for (i in 1:D) {
  wave1imp_bergm[[i]] <- as.matrix.network(bergmImputation$impNets[[i]])
  bergmImputation_numImpEdges[i] <- table(wave1imp_bergm[[i]])[2] #save only number of 1s
}

bergmImputation_numImpEdges #list of edgenumbers in imputed networks
View(wave1imp_bergm[[1]]) #check imputed network no. 1

#save the impuation results
save(wave1imp_bergm, bergmImputation, bergmImputation_numImpEdges, file = "BERGM_Imputation.RData")

#bgof does not work with bergmM (with missings)
#bergm was run again with missings removed (NA = 0)
#The problem with outdegree = 4 is vers similar to the gof results of the ergm model
#for that see temporarily saved > bergm_bgof_test (R script and png file for result), run with main.it = 3000

############## later Waves ###############
load("BERGM_Imputation.RData")

#create lists for the imputed networks of wave 2 and 3
wave2imp_bergm <- list()
wave3imp_bergm <- list()

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
  
  Bekannte_full <- sienaDependent(array(c(wave1imp_bergm[[i]],Bekannte_full_w2_ma),
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
  
  wave2imp_bergm[[i]] <- getNet(Bekannte_full_w2_ma,sims)
  
  # impute wave 3
  
  Bekannte_full <- sienaDependent(array( c(wave2imp_bergm[[i]], Bekannte_full_w3_ma),
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
  
  wave3imp_bergm[[i]] = getNet(Bekannte_full_w3_ma,sims)
  save.image('mi_bergm.RData') 
}

#check results and rename to identify results later:

period1saom_bergm <- period1saom
period2saom_bergm <- period2saom

summary(period1saom_bergm) #cov ok
summary(period2saom_bergm) # cov ok

#check for autocorrelation
period1saom_bergm$ac
period2saom_bergm$ac

# 1: Goodness of Fit
bergm.saom.model.results.gof.p1 <- list()
bergm.saom.model.results.gof.p1[[1]] <- plot(sienaGOF(period1saom_bergm, varName="Bekannte_full", OutdegreeDistribution))
bergm.saom.model.results.gof.p1[[2]] <- plot(sienaGOF(period1saom_bergm, varName="Bekannte_full", IndegreeDistribution))
#View results
bergm.saom.model.results.gof.p1[[1]]
bergm.saom.model.results.gof.p1[[2]]

bergm.saom.model.results.gof.p2 <- list()
bergm.saom.model.results.gof.p2[[1]] <- plot(sienaGOF(period2saom_bergm, varName="Bekannte_full", OutdegreeDistribution))
bergm.saom.model.results.gof.p2[[2]] <- plot(sienaGOF(period2saom_bergm, varName="Bekannte_full", IndegreeDistribution))
#View results
bergm.saom.model.results.gof.p2[[1]]
bergm.saom.model.results.gof.p2[[2]]

save(wave1imp_bergm, wave2imp_bergm, wave3imp_bergm,
     bergmImputation, bergmImputation_numImpEdges, file = "BERGM_Imputation.RData")

########################################################################################################################
####################################### Estimating the analysis models #################################################
########################################################################################################################

############### SAOM Analysis for the Stationary BERGM imputation ##################

#run our model on the D complete data sets in a loop
#output is saved in a list

saomResults_bergm <- list()

for (i in 1:D) {
  cat('Imputation',i,'\n')
  
  Bekannte_full <- sienaDependent(array(c(wave1imp_bergm[[i]], wave2imp_bergm[[i]],
                                          wave3imp_bergm[[i]]), dim = c( 42, 42, 3)))
  
  Data  <- sienaDataCreate(Bekannte_full, Statistik)
  effectsData <- getEffects(Data)
  effectsData <- includeEffects(effectsData,
                                gwespFF, recip, density, outAct)
  effectsData <- includeEffects(effectsData,  sameX,
                                interaction1 =  "Statistik")
  
  estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 312,
                                             n3 = 1000, maxlike = FALSE, cond = FALSE,
                                             lessMem = TRUE) #cond = FALSE wasn't used in Krause (2019) but otherwise the rate parameter is separated in the output
  if (i == 1) {
    saomResults_bergm[[i]] <- siena07ToConvergence(alg = estimation.options,
                                                   dat = Data, eff = effectsData,
                                                   threshold = 0.25)#,
    #nodes = 2)
  } else {
    saomResults_bergm[[i]] <- siena07ToConvergence(alg = estimation.options,
                                                   dat = Data, eff = effectsData,
                                                   threshold = 0.25,
                                                   ans0 = saomResults_bergm[[i - 1]])#, 
    # nodes = 2)
  }
  save.image('mi_bergm.RData') 
}

#check rsults of SAOM for the fisrt imputation
summary(saomResults_bergm[[1]])

############# Combining the Results by Rubin's rules ################

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
  MIResults[,i * 2 - 1] <- saomResults_bergm[[i]]$theta # in original this was $theta, which included the rate parameter, I included it manually above
  MIResults[,i * 2] <-  sqrt(diag(saomResults_bergm[[i]]$covtheta))
}

#does not work if saomResults does not show rate parameter in list 

#Now we get the average covariance structure between the parameters

WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + saomResults_bergm[[i]]$covtheta
}

WDMIs <- (1/D) * WDMIs

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure

finalResults <- as.data.frame(matrix(,npar,2))
names(finalResults) <- c("combinedEstimate", "combinedSE")
rownames(finalResults) <- effectsData$effectName[effectsData$include]
finalResults$combinedEstimate <- rowMeans(MIResults[,seq(1,2*D,2)])
finalResults$combinedSE <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                                  rowVar(MIResults[,seq(1,2*D,2)]))

#rename for identification
finalResults_bergm <- finalResults

#Odds Ratios ergänzen
finalResults_bergm$OddsRatio <- exp(finalResults_bergm$combinedEstimate)


save(wave1imp_bergm, wave2imp_bergm, wave3imp_bergm,
     bergmImputation, bergmImputation_numImpEdges, saomResults_bergm,
     finalResults_bergm_rounded, finalResults_bergm, file = "BERGM_Imputation_Analysis_Results.RData")

#rounding
finalResults_bergm_rounded <- as.data.frame(round(finalResults_bergm_rounded, 3))
View(finalResults_bergm_rounded) #OR und significance fuer rate constants loeschen

#get LaTeX read table
stargazer(finalResults_bergm_rounded, summary = FALSE, title="SAOM Analyse, 1. Welle imputiert mit BERGM",
          rownames = TRUE)

#Significance levels
finalResults_bergm_rounded$sig <- NA
finalResults_bergm_rounded$sig[(abs(finalResults_bergm_rounded$combinedEstimate) - abs(finalResults_bergm_rounded$combinedSE)*1.96)>0]<- "*"
finalResults_bergm_rounded$sig[(abs(finalResults_bergm_rounded$combinedEstimate) - abs(finalResults_bergm_rounded$combinedSE)*2.58)>0]<- "**"
finalResults_bergm_rounded$sig[(abs(finalResults_bergm_rounded$combinedEstimate) - abs(finalResults_bergm_rounded$combinedSE)*3.29)>0]<- "***"

############## Analysis: TERGM Estimation for Stationary SAOM imputation #############

#The analysis is rerun using a separable temporal exponential random graph model. STERGM models formation and dissolution effects separately, hece "separable TERGM".
#here dissolution effects are only included to improve the fit.
#As a model for all data at once did not converge, the model is estimated for both periods separately.

#load the imputed data sets if they are not already in the environment

load("BERGM_Imputation.RData")

#number of imputations (skip if initialized above)
D <- 50

#Convert the SAOM simulations to network for tergm to use
wave1imp_net_bergm <- list()
wave2imp_net_bergm <- list()
wave3imp_net_bergm <- list()

for (i in 1:D) {
  wave1imp_net_bergm[[i]] <- as.network(wave1imp_bergm[[i]])
  wave2imp_net_bergm[[i]] <- as.network(wave2imp_bergm[[i]])
  wave3imp_net_bergm[[i]] <- as.network(wave3imp_bergm[[i]])
}

#Add Statistik attribute to all the simulated networks
set.vertex.attribute(Bekannte_full_net_w1, 'Statistik', Attributes_K1$Statistik_Zero)

for (i in 1:D) {
  set.vertex.attribute(wave1imp_net_bergm[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
  set.vertex.attribute(wave2imp_net_bergm[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
  set.vertex.attribute(wave3imp_net_bergm[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
}


#period 1: wave 1 and 2
Result_tergm_bergm_w12 <- list()

set.seed(234)
for (i in 1:D) {
  cat('Imputation',i,'\n')
  Result_tergm_bergm_w12[[i]] <- tergm(list(wave1imp_net_bergm[[i]], wave2imp_net_bergm[[i]])~
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
summary(Result_tergm_bergm_w12[[1]])


#period 2: wave 2 and 3
Result_tergm_bergm_w23 <- list()

set.seed(234)
for (i in 1:D) {
  cat('Imputation',i,'\n')
  Result_tergm_bergm_w23[[i]] <- tergm(list(wave2imp_net_bergm[[i]], wave3imp_net_bergm[[i]])~
                                         Form(~ edges
                                              + mutual
                                              + gwesp(decay = log(2), fixed = TRUE)
                                              + odegree1.5() #no convergence with odegree1.5()
                                              + nodematch('Statistik'))+
                                         Diss(~ edges
                                              + mutual
                                              + nodematch('Statistik')),
                                       estimate="CMLE")
}
#see result with 1st imputation
summary(Result_tergm_bergm_w23[[1]])

####### Make data set for pooling: wave 1 and 2 ######
npar <- 8 #number of parameters included in the TERGM
TERGM_Results_bergm_combined_w12 <- as.data.frame(matrix(,npar,(4 * D))) #data frame for results

for (i in 1:D) {
  names(TERGM_Results_bergm_combined_w12)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(TERGM_Results_bergm_combined_w12)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  TERGM_Results_bergm_combined_w12[,i * 2 - 1] <- summary(Result_tergm_bergm_w12[[i]])$coefs$"Estimate"
  TERGM_Results_bergm_combined_w12[,i * 2] <-  summary(Result_tergm_bergm_w12[[i]])$coefs$"Std. Error"
}

#Now we get the average covariance structure between the parameters
WDMIs_w12 <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs_w12 <- WDMIs_w12 + Result_tergm_bergm_w12[[i]]$est.cov
}

WDMIs_w12 <- (1/D) * WDMIs_w12

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

finalResults_TERGM_bergm_w12 <- as.data.frame(matrix(,npar,2))
names(finalResults_TERGM_bergm_w12) <- c("combinedEstimate", "combinedSE")
rownames(finalResults_TERGM_bergm_w12) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                            "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_bergm_w12$combinedEstimate <- rowMeans(TERGM_Results_bergm_combined_w12[,seq(1,2*D,2)])
finalResults_TERGM_bergm_w12$combinedSE <- sqrt(diag(WDMIs_w12) + ((D + 1)/D) *
                                                  rowVar(TERGM_Results_bergm_combined_w12[,seq(1,2*D,2)]))

#Odds Ratio hinzufuegen
finalResults_TERGM_bergm_w12$OddsRatio <- exp(finalResults_TERGM_bergm_w12$combinedEstimate)

#rounding
kable(round(finalResults_TERGM_bergm_w12, 3))
finalResults_TERGM_bergm_rounded_w12 <- as.data.frame(round(finalResults_TERGM_bergm_w12, 3))
View(finalResults_TERGM_bergm_rounded_w12)

####### Make data set for pooling wave 2 and 3 ######
npar <- 8 #number of parameters included in the TERGM
TERGM_Results_bergm_combined_w23 <- as.data.frame(matrix(,npar,(4 * D))) #data frame for results

for (i in 1:D) {
  names(TERGM_Results_bergm_combined_w23)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(TERGM_Results_bergm_combined_w23)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  TERGM_Results_bergm_combined_w23[,i * 2 - 1] <- summary(Result_tergm_bergm_w23[[i]])$coefs$"Estimate"
  TERGM_Results_bergm_combined_w23[,i * 2] <-  summary(Result_tergm_bergm_w23[[i]])$coefs$"Std. Error"
}

#Now we get the average covariance structure between the parameters
WDMIs_w23 <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs_w23 <- WDMIs_w23 + Result_tergm_bergm_w23[[i]]$est.cov
}

WDMIs_w23 <- (1/D) * WDMIs_w23

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure
finalResults_TERGM_bergm_w23 <- as.data.frame(matrix(,npar,2))
names(finalResults_TERGM_bergm_w23) <- c("combinedEstimate", "combinedSE")
rownames(finalResults_TERGM_bergm_w23) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.69", "Form~odegree1.5", "Form~nodematch.Statistik",
                                            "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_bergm_w23$combinedEstimate <- rowMeans(TERGM_Results_bergm_combined_w23[,seq(1,2*D,2)])
finalResults_TERGM_bergm_w23$combinedSE <- sqrt(diag(WDMIs_w23) + ((D + 1)/D) *
                                                  rowVar(TERGM_Results_bergm_combined_w23[,seq(1,2*D,2)]))

#Odds Ratios
finalResults_TERGM_bergm_w23$OddsRatio <- exp(finalResults_TERGM_bergm_w23$combinedEstimate)

#Significance
finalResults_TERGM_bergm_w23$OddsRatio_Lower95 <- exp((finalResults_TERGM_bergm_w23$combinedEstimate)-1.96*(finalResults_TERGM_bergm_w23$combinedSE))
finalResults_TERGM_bergm_w23$OddsRatio_Upper95 <- exp((finalResults_TERGM_bergm_w23$combinedEstimate)+1.96*(finalResults_TERGM_bergm_w23$combinedSE))
finalResults_TERGM_bergm_w23$OddsRatio_Lower99 <- exp((finalResults_TERGM_bergm_w23$combinedEstimate)-2.58*(finalResults_TERGM_bergm_w23$combinedSE))
finalResults_TERGM_bergm_w23$OddsRatio_Upper99 <- exp((finalResults_TERGM_bergm_w23$combinedEstimate)+2.58*(finalResults_TERGM_bergm_w23$combinedSE))
finalResults_TERGM_bergm_w23$OddsRatio_Lower99.9 <- exp((finalResults_TERGM_bergm_w23$combinedEstimate)-3.29*(finalResults_TERGM_bergm_w23$combinedSE))
finalResults_TERGM_bergm_w23$OddsRatio_Upper99.9 <- exp((finalResults_TERGM_bergm_w23$combinedEstimate)+3.29*(finalResults_TERGM_bergm_w23$combinedSE))
finalResults_TERGM_bergm_w23$sig <- "-"
finalResults_TERGM_bergm_w23$sig[(finalResults_TERGM_bergm_w23$OddsRatio_Lower95 > 1 & finalResults_TERGM_bergm_w23$OddsRatio_Upper95 >1) |
                                   (finalResults_TERGM_bergm_w23$OddsRatio_Lower95 < 1 & finalResults_TERGM_bergm_w23$OddsRatio_Upper95 <1)] <- "*"
finalResults_TERGM_bergm_w23$sig[(finalResults_TERGM_bergm_w23$OddsRatio_Lower99 > 1 & finalResults_TERGM_bergm_w23$OddsRatio_Upper99 >1) |
                                   (finalResults_TERGM_bergm_w23$OddsRatio_Lower99 < 1 & finalResults_TERGM_bergm_w23$OddsRatio_Upper99 <1)] <- "**"
finalResults_TERGM_bergm_w23$sig[(finalResults_TERGM_bergm_w23$OddsRatio_Lower99.9 > 1 & finalResults_TERGM_bergm_w23$OddsRatio_Upper99.9 >1) |
                                   (finalResults_TERGM_bergm_w23$OddsRatio_Lower99.9 < 1 & finalResults_TERGM_bergm_w23$OddsRatio_Upper99.9 <1)] <- "***"
finalResults_TERGM_bergm_w23 <- finalResults_TERGM_bergm_w23 %>%
  select(-(starts_with("OddsRatio_")))

#rounding
finalResults_TERGM_bergm_rounded_w23 <- as.data.frame(round(finalResults_TERGM_bergm_w23[,1:3], 3))
finalResults_TERGM_bergm_rounded_w23$sig <- finalResults_TERGM_bergm_w23$sig
View(finalResults_TERGM_bergm_rounded_w23)

#table LaTeX ready
library(stargazer)
rownames(finalResults_TERGM_bergm_rounded_w12) <- c("Form edges", "Form mutual", "Form gwesp.fixed.0.69", "Form odegree1.5", "Form nodematch.Statistik",
                                                    "Diss edges", "Diss mutual", "Diss nodematch.Statistik")
stargazer(finalResults_TERGM_bergm_rounded_w12, summary = FALSE, title="TERGM Anayse, 1. Welle imputiert mit BERGM, Periode 1",
          rownames = TRUE)

#combined
finalResults_TERGM_bergm_rounded_w12_and_w23 <- finalResults_TERGM_bergm_rounded_w12
finalResults_TERGM_bergm_rounded_w12_and_w23$OddsRatio <- round(exp(finalResults_TERGM_bergm_rounded_w12$combinedEstimate), 3)
finalResults_TERGM_bergm_rounded_w12_and_w23$combinedEstimate23 <- finalResults_TERGM_bergm_rounded_w23$combinedEstimate
finalResults_TERGM_bergm_rounded_w12_and_w23$combinedSE23 <- finalResults_TERGM_bergm_rounded_w23$combinedSE
finalResults_TERGM_bergm_rounded_w12_and_w23$OddsRatio23 <- round(exp(finalResults_TERGM_bergm_rounded_w23$combinedEstimate), 3)

stargazer(finalResults_TERGM_bergm_rounded_w12_and_w23, summary = FALSE, title="TERGM Anayse, 1. Welle imputiert mit BERGM",
          rownames = TRUE, column.labels= c("Periode 1", "Periode 2"), column.separate = c(3,3))

save(wave1imp_net_bergm, wave2imp_net_bergm, wave3imp_net_bergm,
     finalResults_TERGM_bergm_w12, finalResults_TERGM_bergm_rounded_w12, TERGM_Results_bergm_combined_w12,
     Result_tergm_bergm_w12, finalResults_TERGM_bergm_w23, finalResults_TERGM_bergm_rounded_w23,
     TERGM_Results_bergm_combined_w23, Result_tergm_bergm_w23,
     file = "BERGM_TERGM_Analysis_Results.RData")

################### Assessing goodness of fit #############################
gof_tergm_bergm_w12 <- gof(Result_tergm_bergm_w12[[6]])
gof_tergm_bergm_w12 
plot(gof_tergm_bergm_w12)

gof_tergm_bergm_w23 <- gof(Result_tergm_bergm_w23[[6]])
gof_tergm_bergm_w23 
plot(gof_tergm_bergm_w23)

mcmc.diagnostics(Result_tergm_bergm_w12[[6]])

#Significance

finalResults_TERGM_bergm_rounded_w12_and_w23$sig_w12 <- NA
finalResults_TERGM_bergm_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedSE)*1.96)>0]<- "*"
finalResults_TERGM_bergm_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedSE)*2.58)>0]<- "**"
finalResults_TERGM_bergm_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedSE)*3.29)>0]<- "***"

finalResults_TERGM_bergm_rounded_w12_and_w23$sig_w23 <- NA
finalResults_TERGM_bergm_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedSE23)*1.96)>0]<- "*"
finalResults_TERGM_bergm_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedSE23)*2.58)>0]<- "**"
finalResults_TERGM_bergm_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_bergm_rounded_w12_and_w23$combinedSE23)*3.29)>0]<- "***"



