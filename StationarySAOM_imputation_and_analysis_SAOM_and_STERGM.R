#Code for the imputation of HeiNet with Stationary SAOM and for analysis with SAOM and with STERGM

#An equivalent was applied for the BERGM, ERGM and SAOM internal imputation.
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

#get the data
load("mydata.RData")
#Bekannte_full_w1_ma, Bekannte_full_w2_ma and Bekannte_full_w3_ma are adjecency matrices in the matrix format.
#They have a size of 42x42 and represent the the waves of relational data collected, missings are NA.
#The data also includes attribute data about the tutorial group in statistics each individual participated in in the data.table Attributes_K1.
#Missings in the attributes are not imputed.

#DEfine Function that will hopefully ensure convergence in RSiena
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
wave1imp <- list() #create list to store the imputed network
D <- 50 #number of imputations


## Imputation with a stationary SAOM
Bekannte_full <- sienaDependent(array(c(Bekannte_full_w1_ma, Bekannte_full_w1_ma), dim = c( 42, 42, 2)) ,
                                allowOnly = FALSE)

Statistik_v <- as.vector(Attributes_K1$Statistik) #Attribute Data

Statistik  <- coCovar(as.vector(Statistik_v),centered=FALSE)
w2 <- coDyadCovar(Bekannte_full_w2_ma) #using the second wave as a dyadic covariate for the imputation of wave 1

Data.stationary <- sienaDataCreate(Bekannte_full,Statistik,w2)
effects.stationary <- getEffects(Data.stationary)
effects.stationary <- includeEffects(effects.stationary, recip, density, gwespFF, outAct)

effects.stationary <- includeEffects(effects.stationary, sameX,
                                     interaction1 = "Statistik")
effects.stationary <- includeEffects(effects.stationary, X, name = "Bekannte_full",
                                     interaction1 = "w2")

effects.stationary <- setEffect(effects.stationary, Rate, initialValue = 50,
                                name = "Bekannte_full", fix = TRUE, type = "rate")

#Now we can estimate the stationary SAOM with Methods of Moments (MoM) estimation. One converged estimate will suffice.

estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE,
                                           seed = 2214,
                                           n3 = 1000, maxlike = FALSE,
                                           cond = FALSE) #n3 =number of datasets simulated

period0saom <- siena07ToConvergence(alg = estimation.options,
                                    dat = Data.stationary,
                                    eff = effects.stationary, threshold = 0.25)
#The model converged in the second period

#Now change the RSiena algorithm to impuatation

imputation.options <- sienaAlgorithmCreate(seed = 13848,
                                           useStdInits = FALSE, 
                                           maxlike = TRUE,
                                           cond = FALSE, 
                                           nsub = 0,
                                           simOnly = TRUE,
                                           n3 = 10)

#Goodness of fit should ideally be evaluated before using the model for imputation (see below).
#In this case it is compromised for comparability with the other imputation mechanisms and with the STERGM analysis method.

#Show results of the above model
period0saom
#t-ratio of the basic rate parameter is so high and SE is NA, because the effect is fixed

#check for autocorrelation
period0saom$ac

#check for collinearity
summary(period0saom)
#gwespFF und gwespFB scheinen eine relativ hohe Kolinearität zu haben

# 1: Goodness of Fit
saom.model.results.gof <- list()
saom.model.results.gof[[1]] <- plot(sienaGOF(period0saom, varName="Bekannte_full", OutdegreeDistribution))
saom.model.results.gof[[2]] <- plot(sienaGOF(period0saom, varName="Bekannte_full", IndegreeDistribution))
#View results
saom.model.results.gof[[1]]
saom.model.results.gof[[2]]
#ggsave(file="temporarily saved/GoF_1_Selection.pdf",width=50,height=35,units="cm")


########################################################################################################################

### To prepare the imputation two steps are necessary
#1) One tie needs to be different between the waves (wave one is treated as the starting and the end wave for the imputation). It will be correctly changed back by the mL algorithm
#2) All observed ties (except the changed one) will be set to structurally fixed values (1 = 11; 0 = 10)

set.seed(142)
for (i in 1:D) {
  cat('imputation',i,'\n')
  
  n1 <- Bekannte_full_w1_ma
  n1 <- n1 + 10
  diag(n1) <- 0
  n2 <- n1
  tieList <- c(1:(nrow(n1)**2))[c(n1 == 11)]
  tieList <- tieList[!is.na(tieList)]
  
  changedTie <- sample(tieList,1)
  
  n1[changedTie] <- 0
  n2[changedTie] <- 1
  
  Bekannte_full <- sienaDependent(array(c(n1,n2), dim = c( 42, 42, 2)),
                                  allowOnly = FALSE )
  
  Data.stationary <- sienaDataCreate(Bekannte_full,Statistik,w2)
  
  sims <- siena07(imputation.options, data = Data.stationary,
                  effects = effects.stationary,
                  prevAns = period0saom,
                  returnDeps = TRUE)$sims[[10]][[1]]
  
  wave1imp[[i]] = getNet(Bekannte_full_w1_ma,sims)
  
}

### Worked!!

View(wave1imp[[1]]) #test

############## later Waves ###############

#create lists for the imputed networks of wave 2 and 3
wave2imp <- list()
wave3imp <- list()

#For each wave estimation is done WIth Method of Moments (MoM) and imputation with Maximum Likelihood (ML)
#Imputation is repeated as often as the first wave had been imputed (D times)

set.seed(1307)
for (i in 1:D) {
  
  cat('imputation',i,'\n')
  
  #impute wave2 with the SAOM internal mechanism
  
  Bekannte_full <- sienaDependent(array(c(wave1imp[[i]],Bekannte_full_w2_ma),
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
  
  wave2imp[[i]] <- getNet(Bekannte_full_w2_ma,sims)
  
  # impute wave 3
  
  Bekannte_full <- sienaDependent(array( c(wave2imp[[i]], Bekannte_full_w3_ma),
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
  
  wave3imp[[i]] = getNet(Bekannte_full_w3_ma,sims)
  save.image('mi.RData') 
}

#check results:

period1saom
period2saom

summary(period1saom)
summary(period2saom)

#check for autocorrelation
period1saom$ac
period2saom$ac

#check for collinearity
summary(period0saom) 
summary(period1saom)
summary(period2saom) 

# 1: Goodness of Fit
saom.model.results.gof.p1 <- list()
saom.model.results.gof.p1[[1]] <- plot(sienaGOF(period1saom, varName="Bekannte_full", OutdegreeDistribution))
saom.model.results.gof.p1[[2]] <- plot(sienaGOF(period1saom, varName="Bekannte_full", IndegreeDistribution))
#View results
saom.model.results.gof.p1[[1]]
saom.model.results.gof.p1[[2]]

saom.model.results.gof.p2 <- list()
saom.model.results.gof.p2[[1]] <- plot(sienaGOF(period2saom, varName="Bekannte_full", OutdegreeDistribution))
saom.model.results.gof.p2[[2]] <- plot(sienaGOF(period2saom, varName="Bekannte_full", IndegreeDistribution))
#View results
saom.model.results.gof.p2[[1]]
saom.model.results.gof.p2[[2]]

#To call certain Imputations: For Example 4th imputation of wave 2:
View(wave2imp[[4]])

######## Renaming the objects to be exported ###########
wave1imp_StatSAOM <- wave1imp
wave2imp_StatSAOM <- wave2imp
wave3imp_StatSAOM <- wave3imp

stat.saom.model.results.gof.p1 <- saom.model.results.gof.p1
stat.saom.model.results.gof.p2 <- saom.model.results.gof.p2
stat.saom.model.results.gof.p0 <- saom.model.results.gof

period0StatSAOM <- period0saom
period1StatSAOM <- period1saom
period2StatSAOM <- period2saom

#save the model results
save(wave1imp_StatSAOM, wave2imp_StatSAOM, wave3imp_StatSAOM,
     stat.saom.model.results.gof.p1, stat.saom.model.results.gof.p2, stat.saom.model.results.gof.p0,
     period0StatSAOM, period1StatSAOM, period2StatSAOM, file = "StatSAOM_Imputation.RData")


########################################################################################################################
####################################### Estimating the analysis models #################################################
########################################################################################################################

############### SAOM Analysis for the Stationary SAOM imputation ##################

load("StatSAOM_Imputation.RData")

#run model on the D complete data sets in a loop
#output is saved in a list

saomResults <- list()

for (i in 1:D) {
  cat('Imputation',i,'\n')
  
  Bekannte_full <- sienaDependent(array(c(wave1imp[[i]], wave2imp[[i]],
                                          wave3imp[[i]]), dim = c( 42, 42, 3)))
  
  Data  <- sienaDataCreate(Bekannte_full, Statistik)
  effectsData <- getEffects(Data)
  effectsData <- includeEffects(effectsData,
                                gwespFF, recip, density, outAct)
  effectsData <- includeEffects(effectsData,  sameX,
                                interaction1 =  "Statistik")
  
  estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 312,
                                             n3 = 1000, maxlike = FALSE, cond = FALSE,
                                             lessMem = FALSE) #cond = FALSE wasn't used in Krause (2019) but otherwise the rate parameter is separated in the output
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
  save.image('mi_StatSAOM.RData') 
}

saomResults[[1]] #estimation using the first imputation: Range from 1-50 (if D = 50)

#colliearity check
summary(saomResults[[1]])

#check for autocorrelation
saomResults[[1]]$ac

#gof
saom.results.gof.resultStatSAOM <- list()
saom.results.gof.resultStatSAOM[[1]] <- plot(sienaGOF(saomResults[[1]], varName="Bekannte_full", OutdegreeDistribution))
saom.results.gof.resultStatSAOM[[2]] <- plot(sienaGOF(saomResults[[1]], varName="Bekannte_full", IndegreeDistribution))
saom.results.gof.resultStatSAOM[[1]]
saom.results.gof.resultStatSAOM[[2]]

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
  MIResults[,i * 2 - 1] <- saomResults[[i]]$theta # in original this was $theta, which included the rate parameter, I included it manually above
  MIResults[,i * 2] <-  sqrt(diag(saomResults[[i]]$covtheta))
}

#does not work if saomResults does not show rate parameter in list 

#Now we get the average covariance structure between the parameters

WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + saomResults[[i]]$covtheta
}

WDMIs <- (1/D) * WDMIs

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure

finalResultsStatSAOM <- as.data.frame(matrix(,npar,2))
names(finalResultsStatSAOM) <- c("combinedEstimate", "combinedSE")
rownames(finalResultsStatSAOM) <- effectsData$effectName[effectsData$include]
finalResultsStatSAOM$combinedEstimate <- rowMeans(MIResults[,seq(1,2*D,2)])
finalResultsStatSAOM$combinedSE <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                                          rowVar(MIResults[,seq(1,2*D,2)]))

#Odds Ratios
finalResultsStatSAOM$OddsRatio <- exp(finalResultsStatSAOM$combinedEstimate)

#rounding
finalResultsStatSAOM_rounded <- as.data.frame(round(finalResultsStatSAOM, 3))
View(finalResultsStatSAOM_rounded)
save.image('mi.RData')

MIResults_StatSAOM <- MIResults

#final significance
significance_StatSAOM_Results <- ifelse((abs(finalResultsStatSAOM_rounded$combinedEstimate/finalResultsStatSAOM_rounded$combinedSE)) > 1.96, TRUE, FALSE) #MIND: Absolute value
significance_StatSAOM_Results #all significant

#save the model results
save(saomResults, finalResultsStatSAOM, effectsData, MIResults_StatSAOM, npar, WDMIs,
     finalResultsStatSAOM, finalResultsStatSAOM_rounded, file = "StatSAOM_Imputation_Analysis_Results.RData")

#Odds Ratios
finalResultsStatSAOM_rounded$OddsRatio <- exp(finalResultsStatSAOM_rounded$combinedEstimate)

#rounding
finalResultsStatSAOM_rounded <- as.data.frame(round(finalResultsStatSAOM_rounded, 3))
View(finalResultsStatSAOM_rounded) #OR und significance fuer rate constants loeschen


#getting a LaTeX ready results table
stargazer(finalResultsStatSAOM_rounded, summary = FALSE, title="SAOM Analyse, 1. Welle imputiert mit Stat.SAOM",
          rownames = TRUE)

#Signicance levels
finalResultsStatSAOM_rounded$sig <- NA
finalResultsStatSAOM_rounded$sig[(abs(finalResultsStatSAOM_rounded$combinedEstimate) - abs(finalResultsStatSAOM_rounded$combinedSE)*1.96)>0]<- "*"
finalResultsStatSAOM_rounded$sig[(abs(finalResultsStatSAOM_rounded$combinedEstimate) - abs(finalResultsStatSAOM_rounded$combinedSE)*2.58)>0]<- "**"
finalResultsStatSAOM_rounded$sig[(abs(finalResultsStatSAOM_rounded$combinedEstimate) - abs(finalResultsStatSAOM_rounded$combinedSE)*3.29)>0]<- "***"


############## Analysis: TERGM Estimation for Stationary SAOM imputation #############

#The analysis is rerun using a separable temporal exponential random graph model. STERGM models formation and dissolution effects separately, hece "separable TERGM".
#here dissolution effects are only included to improve the fit.
#As a model for all data at once did not converge, the model is estimated for both periods separately.

#load the imputed data sets if they are not already in the environment
load("StatSAOM_Imputation.RData")
load("StatSAOM_Imputation_Analysis_Results.RData")

#number of imputations (skip if initialized above)
D <- 50

#Convert the SAOM simulations to network for tergm to use (statnet framework)
wave1imp_net_SAOMStat <- list()
wave2imp_net_SAOMStat <- list()
wave3imp_net_SAOMStat <- list()

for (i in 1:D) {
  wave1imp_net_SAOMStat[[i]] <- as.network(wave1imp_StatSAOM[[i]]) 
  wave2imp_net_SAOMStat[[i]] <- as.network(wave2imp_StatSAOM[[i]])
  wave3imp_net_SAOMStat[[i]] <- as.network(wave3imp_StatSAOM[[i]])
}

#Add "Statistik" Attribute to all the simulated networks
set.vertex.attribute(Bekannte_full_net_w1, 'Statistik', Attributes_K1$Statistik_Zero)

for (i in 1:D) {
  set.vertex.attribute(wave1imp_net_SAOMStat[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
  set.vertex.attribute(wave2imp_net_SAOMStat[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
  set.vertex.attribute(wave3imp_net_SAOMStat[[i]], 'Statistik', Attributes_K1$Statistik_Zero)
}


#period 1: wave 1 und 2
Result_tergm_SAOMStat_w12 <- list()

for (i in 1:D) {
  cat('Imputation',i,'\n')
  Result_tergm_SAOMStat_w12[[i]] <- tergm(list(wave1imp_net_SAOMStat[[i]], wave2imp_net_SAOMStat[[i]])~
                                            Form(~ edges
                                                 + mutual
                                                 + gwesp(decay = log(2), fixed = TRUE)
                                                 + odegree1.5() 
                                                 + nodematch('Statistik'))+
                                            Diss(~ edges
                                                 + mutual
                                                 + nodematch('Statistik')),
                                          estimate="CMLE")
  #control = control.tergm(CMLE.MCMC.interval = 4096)) #4096 is 4 times 1024 (the default)
}
#see result with 1st imputation
summary(Result_tergm_SAOMStat_w12[[1]])
#accessing specific results of other imputations e.g., no 6
Result_tergm_SAOMStat_w12[[6]]$sample

#period 2: wave 2 und 3
Result_tergm_SAOMStat_w23 <- list()

for (i in 1:D) {
  cat('Imputation',i,'\n')
  Result_tergm_SAOMStat_w23[[i]] <- tergm(list(wave2imp_net_SAOMStat[[i]], wave3imp_net_SAOMStat[[i]])~
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
summary(Result_tergm_SAOMStat_w23[[1]])

#display the structure of an object
str(summary(Result_tergm_SAOMStat_w12[[1]]))

#How to extract the individual components of the results? Here e.g., imputation 4
summary(Result_tergm_SAOMStat_w12[[4]])$coefs #all coefs
summary(Result_tergm_SAOMStat_w12[[4]])$coefs$"Std. Error" #only extract SE
summary(Result_tergm_SAOMStat_w12[[4]])$coefs$"Pr(>|z|)"
summary(Result_tergm_SAOMStat_w12[[4]])$coefs$"Estimate"

Result_tergm_SAOMStat_w12[[4]]$est.cov #covariance matrix of the estimates

####### Make data set for pooling wave 1 and 2 ######
npar <- 8 #number of parameters included in the TERGM
TERGM_Results_StatSAOM_combined_w12 <- as.data.frame(matrix(,npar,(4 * D))) #data frame for results

for (i in 1:D) {
  names(TERGM_Results_StatSAOM_combined_w12)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(TERGM_Results_StatSAOM_combined_w12)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  TERGM_Results_StatSAOM_combined_w12[,i * 2 - 1] <- summary(Result_tergm_SAOMStat_w12[[i]])$coefs$"Estimate"
  TERGM_Results_StatSAOM_combined_w12[,i * 2] <-  summary(Result_tergm_SAOMStat_w12[[i]])$coefs$"Std. Error"
}

#Now we get the average covariance structure between the parameters
WDMIs_w12 <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs_w12 <- WDMIs_w12 + Result_tergm_SAOMStat_w12[[i]]$est.cov
}

WDMIs_w12 <- (1/D) * WDMIs_w12

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

finalResults_TERGM_StatSAOM_w12 <- as.data.frame(matrix(,npar,2))
names(finalResults_TERGM_StatSAOM_w12) <- c("combinedEstimate", "combinedSE")
rownames(finalResults_TERGM_StatSAOM_w12) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                               "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_StatSAOM_w12$combinedEstimate <- rowMeans(TERGM_Results_StatSAOM_combined_w12[,seq(1,2*D,2)])
finalResults_TERGM_StatSAOM_w12$combinedSE <- sqrt(diag(WDMIs_w12) + ((D + 1)/D) *
                                                     rowVar(TERGM_Results_StatSAOM_combined_w12[,seq(1,2*D,2)]))

#Odds Ratios
finalResults_TERGM_StatSAOM_w12$OddsRatio <- exp(finalResults_TERGM_StatSAOM_w12$combinedEstimate)

#rounding
finalResults_TERGM_StatSAOM_rounded_w12 <- as.data.frame(round(finalResults_TERGM_StatSAOM_w12, 3))
View(finalResults_TERGM_StatSAOM_rounded_w12)

####### Make data set for pooling wave 2 and 3 ######
npar <- 8 #number of parameters included in the TERGM
TERGM_Results_StatSAOM_combined_w23 <- as.data.frame(matrix(,npar,(4 * D))) #data frame for results

for (i in 1:D) {
  names(TERGM_Results_StatSAOM_combined_w23)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(TERGM_Results_StatSAOM_combined_w23)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  TERGM_Results_StatSAOM_combined_w23[,i * 2 - 1] <- summary(Result_tergm_SAOMStat_w23[[i]])$coefs$"Estimate"
  TERGM_Results_StatSAOM_combined_w23[,i * 2] <-  summary(Result_tergm_SAOMStat_w23[[i]])$coefs$"Std. Error"
}

#Now we get the average covariance structure between the parameters
WDMIs_w23 <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs_w23 <- WDMIs_w23 + Result_tergm_SAOMStat_w23[[i]]$est.cov
}

WDMIs_w23 <- (1/D) * WDMIs_w23

#Using Rubin’s Rules we combine the parameters and standard errors and complete the procedure
finalResults_TERGM_StatSAOM_w23 <- as.data.frame(matrix(,npar,2))
names(finalResults_TERGM_StatSAOM_w23) <- c("combinedEstimate", "combinedSE")
rownames(finalResults_TERGM_StatSAOM_w23) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                               "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_StatSAOM_w23$combinedEstimate <- rowMeans(TERGM_Results_StatSAOM_combined_w23[,seq(1,2*D,2)])
finalResults_TERGM_StatSAOM_w23$combinedSE <- sqrt(diag(WDMIs_w23) + ((D + 1)/D) *
                                                     rowVar(TERGM_Results_StatSAOM_combined_w23[,seq(1,2*D,2)]))


#Odds Ratios
finalResults_TERGM_StatSAOM_w23$OddsRatio <- exp(finalResults_TERGM_StatSAOM_w23$combinedEstimate)

#rounding
finalResults_TERGM_StatSAOM_rounded_w23 <- as.data.frame(round(finalResults_TERGM_StatSAOM_w23, 3))
View(finalResults_TERGM_StatSAOM_rounded_w23)

save(wave1imp_net_SAOMStat, wave2imp_net_SAOMStat, wave3imp_net_SAOMStat,
     finalResults_TERGM_StatSAOM_w12, finalResults_TERGM_StatSAOM_rounded_w12, TERGM_Results_StatSAOM_combined_w12,
     Result_tergm_SAOMStat_w12, finalResults_TERGM_StatSAOM_w23, finalResults_TERGM_StatSAOM_rounded_w23,
     TERGM_Results_StatSAOM_combined_w23, Result_tergm_SAOMStat_w23,
     file = "Stat_SAOM_TERGM_Analysis_Results.RData")

################### Assessing goodness of fit #############################
gof_tergm_SAOMStat_w12 <- gof(Result_tergm_SAOMStat_w12[[6]])
gof_tergm_SAOMStat_w12 
plot(gof_tergm_SAOMStat_w12)

mcmc.diagnostics(Result_tergm_SAOMStat_w12[[6]])
mcmc.diagnostics(Result_tergm_SAOMStat_w12[[16]])

gof_tergm_SAOMStat_w23 <- gof(Result_tergm_SAOMStat_w23[[6]])
gof_tergm_SAOMStat_w23 

#make LaTeX ready table
rownames(finalResults_TERGM_StatSAOM_rounded_w12) <- c("F edges", "F mutual", "F gwesp(0.69)", "F odegree1.5", "F Fokus Statistik",
                                                       "D edges", "D mutual", "D Fokus Statistik")

#combined
finalResults_TERGM_StatSAOM_rounded_w12_and_w23 <- finalResults_TERGM_StatSAOM_rounded_w12
finalResults_TERGM_StatSAOM_rounded_w12_and_w23$OddsRatio <- round(exp(finalResults_TERGM_StatSAOM_rounded_w12$combinedEstimate), 3)
finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedEstimate23 <- finalResults_TERGM_StatSAOM_rounded_w23$combinedEstimate
finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedSE23 <- finalResults_TERGM_StatSAOM_rounded_w23$combinedSE
finalResults_TERGM_StatSAOM_rounded_w12_and_w23$OddsRatio23 <- round(exp(finalResults_TERGM_StatSAOM_rounded_w23$combinedEstimate), 3)

stargazer(finalResults_TERGM_StatSAOM_rounded_w12_and_w23, summary = FALSE, title="TERGM Anayse, 1. Welle imputiert mit Stationary SAOM",
          rownames = TRUE, column.labels= c("Periode 1", "Periode 2"), column.separate = c(3,3))


#Significance levels

finalResults_TERGM_StatSAOM_rounded_w12_and_w23$sig_w12 <- NA
finalResults_TERGM_StatSAOM_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedSE)*1.96)>0]<- "*"
finalResults_TERGM_StatSAOM_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedSE)*2.58)>0]<- "**"
finalResults_TERGM_StatSAOM_rounded_w12_and_w23$sig_w12[(abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedEstimate) - abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedSE)*3.29)>0]<- "***"

finalResults_TERGM_StatSAOM_rounded_w12_and_w23$sig_w23 <- NA
finalResults_TERGM_StatSAOM_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedSE23)*1.96)>0]<- "*"
finalResults_TERGM_StatSAOM_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedSE23)*2.58)>0]<- "**"
finalResults_TERGM_StatSAOM_rounded_w12_and_w23$sig_w23[(abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedEstimate23) - abs(finalResults_TERGM_StatSAOM_rounded_w12_and_w23$combinedSE23)*3.29)>0]<- "***"


