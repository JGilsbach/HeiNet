#analysis with SAOM and with STERGM of HeiNet under complete null imputation (all waves NA = 0)

#for comparison with the multiply imputed analyses

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

#make an NA = 0 version of wave 1
Bekannte_full_w1_ma_null <- Bekannte_full_w1_ma
Bekannte_full_w1_ma_null[is.na(Bekannte_full_w1_ma_null)] <- 0
Bekannte_full_net_w1_null <- as.network(Bekannte_full_w1_ma_null)

#make an NA = 0 version of wave 2
Bekannte_full_w2_ma_null <- Bekannte_full_w2_ma
Bekannte_full_w2_ma_null[is.na(Bekannte_full_w2_ma_null)] <- 0
Bekannte_full_net_w2_null <- as.network(Bekannte_full_w2_ma_null)

#make an NA = 0 version of wave 3
Bekannte_full_w3_ma_null <- Bekannte_full_w3_ma
Bekannte_full_w3_ma_null[is.na(Bekannte_full_w3_ma_null)] <- 0
Bekannte_full_net_w3_null <- as.network(Bekannte_full_w3_ma_null)

#Add Statistik attribute to all the simulated networks
set.vertex.attribute(Bekannte_full_net_w1_null, 'Statistik', Attributes_K1$Statistik_Zero)
set.vertex.attribute(Bekannte_full_net_w2_null, 'Statistik', Attributes_K1$Statistik_Zero)
set.vertex.attribute(Bekannte_full_net_w3_null, 'Statistik', Attributes_K1$Statistik_Zero)

############## Analysis: TERGM Estimation ##############

#period 1: wave 1 and 2
Result_tergm_nullImputation_w12 <- tergm(list(Bekannte_full_net_w1_null, Bekannte_full_net_w2_null)~
                                           Form(~ edges
                                                + mutual
                                                + gwesp(decay = log(2), fixed = TRUE)
                                                + odegree1.5() 
                                                + nodematch('Statistik'))+
                                           Diss(~ edges
                                                + mutual
                                                + nodematch('Statistik')),
                                         estimate="CMLE")

summary(Result_tergm_nullImputation_w12)

Result_tergm_nullImputation_w23 <- tergm(list(Bekannte_full_net_w2_null, Bekannte_full_net_w3_null)~
                                           Form(~ edges
                                                + mutual
                                                + gwesp(decay = log(2), fixed = TRUE)
                                                + odegree1.5() 
                                                + nodematch('Statistik'))+
                                           Diss(~ edges
                                                + mutual
                                                + nodematch('Statistik')),
                                         estimate="CMLE")

summary(Result_tergm_nullImputation_w23)

#goodness of fit and mcmc diagnostics
Result_tergm_nullImputation_w12_gof <- gof(Result_tergm_nullImputation_w12)
Result_tergm_nullImputation_w12_gof 
plot(Result_tergm_nullImputation_w12_gof)

mcmc.diagnostics(Result_tergm_nullImputation_w12)

Result_tergm_nullImputation_w23_gof <- gof(Result_tergm_nullImputation_w23)
Result_tergm_nullImputation_w23_gof 
plot(Result_tergm_nullImputation_w23_gof)

mcmc.diagnostics(Result_tergm_nullImputation_w23)

#save results
save(Result_tergm_nullImputation_w12, Result_tergm_nullImputation_w23,
     Bekannte_full_net_w1_null, Bekannte_full_net_w2_null, Bekannte_full_net_w3_null,
     file = "CompleteNullImp_TERGM_Analysis_Results.RData")

#final dataframe with results
finalResults_TERGM_nullImputation_w12 <- as.data.frame(matrix(,8,2))
names(finalResults_TERGM_nullImputation_w12) <- c("Estimate", "SE")
rownames(finalResults_TERGM_nullImputation_w12) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                                     "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_nullImputation_w12$Estimate <- summary(Result_tergm_nullImputation_w12)$coefficients[1:8]
finalResults_TERGM_nullImputation_w12$SE <- summary(Result_tergm_nullImputation_w12)$asyse

#Odds Ratios
finalResults_TERGM_nullImputation_w12$OddsRatio <- exp(finalResults_TERGM_nullImputation_w12$Estimate)

#rounding
finalResults_TERGM_nullImputation_rounded_w12 <- as.data.frame(round(finalResults_TERGM_nullImputation_w12, 3))
View(finalResults_TERGM_nullImputation_rounded_w12)

#period 2: wave 2 and 3
finalResults_TERGM_nullImputation_w23 <- as.data.frame(matrix(,8,2))
names(finalResults_TERGM_nullImputation_w23) <- c("Estimate", "SE")
rownames(finalResults_TERGM_nullImputation_w23) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                                     "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_nullImputation_w23$Estimate <- summary(Result_tergm_nullImputation_w23)$coefficients[1:8]
finalResults_TERGM_nullImputation_w23$SE <- summary(Result_tergm_nullImputation_w23)$asyse

#Odds Ratios
finalResults_TERGM_nullImputation_w23$OddsRatio <- exp(finalResults_TERGM_nullImputation_w23$Estimate)

#rounding
finalResults_TERGM_nullImputation_rounded_w23 <- as.data.frame(round(finalResults_TERGM_nullImputation_w23, 3))
View(finalResults_TERGM_nullImputation_rounded_w23)

#### together: combine results of both periods
finalResults_TERGM_nullImputation_rounded_w12_w23 <- finalResults_TERGM_nullImputation_rounded_w12
finalResults_TERGM_nullImputation_rounded_w12_w23$Estimate23 <- finalResults_TERGM_nullImputation_rounded_w23$Estimate
finalResults_TERGM_nullImputation_rounded_w12_w23$SE23 <- finalResults_TERGM_nullImputation_rounded_w23$SE
finalResults_TERGM_nullImputation_rounded_w12_w23$OddsRatio23 <- finalResults_TERGM_nullImputation_rounded_w23$OddsRatio
View(finalResults_TERGM_nullImputation_rounded_w12_w23)
finalResults_TERGM_nullImputation_rounded_w12_w23 <- as.data.frame(finalResults_TERGM_nullImputation_rounded_w12_w23)

#get LaTeX ready table
stargazer(finalResults_TERGM_nullImputation_rounded_w12_w23, summary = FALSE, title="TERGM analysis, all waves null imputed",
          rownames = TRUE, column.labels= c("Periode 1", "Periode 2"), column.separate = c(3,3))


#significance levels
finalResults_TERGM_nullImputation_rounded_w12_w23$sig_w12 <- NA
finalResults_TERGM_nullImputation_rounded_w12_w23$sig_w12[(abs(finalResults_TERGM_nullImputation_rounded_w12_w23$Estimate) - abs(finalResults_TERGM_nullImputation_rounded_w12_w23$SE)*1.96)>0]<- "*"
finalResults_TERGM_nullImputation_rounded_w12_w23$sig_w12[(abs(finalResults_TERGM_nullImputation_rounded_w12_w23$Estimate) - abs(finalResults_TERGM_nullImputation_rounded_w12_w23$SE)*2.58)>0]<- "**"
finalResults_TERGM_nullImputation_rounded_w12_w23$sig_w12[(abs(finalResults_TERGM_nullImputation_rounded_w12_w23$Estimate) - abs(finalResults_TERGM_nullImputation_rounded_w12_w23$SE)*3.29)>0]<- "***"

finalResults_TERGM_nullImputation_rounded_w12_w23$sig_w23 <- NA
finalResults_TERGM_nullImputation_rounded_w12_w23$sig_w23[(abs(finalResults_TERGM_nullImputation_rounded_w12_w23$Estimate23) - abs(finalResults_TERGM_nullImputation_rounded_w12_w23$SE23)*1.96)>0]<- "*"
finalResults_TERGM_nullImputation_rounded_w12_w23$sig_w23[(abs(finalResults_TERGM_nullImputation_rounded_w12_w23$Estimate23) - abs(finalResults_TERGM_nullImputation_rounded_w12_w23$SE23)*2.58)>0]<- "**"
finalResults_TERGM_nullImputation_rounded_w12_w23$sig_w23[(abs(finalResults_TERGM_nullImputation_rounded_w12_w23$Estimate23) - abs(finalResults_TERGM_nullImputation_rounded_w12_w23$SE23)*3.29)>0]<- "***"


############ SAOM analysis for null imputation ################

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

#prepare attribute
Statistik_v <- as.vector(Attributes_K1$Statistik)
Statistik  <- coCovar(as.vector(Statistik_v),centered=FALSE)


Bekannte_full <- sienaDependent(array(c(Bekannte_full_w1_ma_null, Bekannte_full_w2_ma_null,
                                        Bekannte_full_w3_ma_null), dim = c( 42, 42, 3)))

Data  <- sienaDataCreate(Bekannte_full, Statistik)
effectsData <- getEffects(Data)
effectsData <- includeEffects(effectsData,
                              gwespFF, recip, density, outAct, outIso) #only converges with additional outIso
effectsData <- includeEffects(effectsData,  sameX,
                              interaction1 =  "Statistik")

estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 312,
                                           n3 = 1000, maxlike = FALSE, cond = FALSE,
                                           lessMem = FALSE)

saomResults_null_imp <- siena07ToConvergence(alg = estimation.options,
                                             dat = Data, eff = effectsData,
                                             threshold = 0.25)

saomResults_null_imp




#colliearity check
summary(saomResults_null_imp)

#check for autocorrelation
saomResults_null_imp$ac

#gof
saom.results.gof.result_Null <- list()
saom.results.gof.result_Null[[1]] <- plot(sienaGOF(saomResults_null_imp, varName="Bekannte_full", OutdegreeDistribution))
saom.results.gof.result_Null[[2]] <- plot(sienaGOF(saomResults_null_imp, varName="Bekannte_full", IndegreeDistribution))
saom.results.gof.result_Null[[1]] 
saom.results.gof.result_Null[[2]] 

#save results
save(saomResults_null_imp,
     Bekannte_full_net_w1_null, Bekannte_full_net_w2_null, Bekannte_full_net_w3_null, saom.results.gof.result_Null,
     file = "CompleteNullImp_SAOM_Analysis_Results.RData")

finalResults_null_imp <- as.data.frame(matrix(,8,2))
names(finalResults_null_imp) <- c("Estimate", "SE")
rownames(finalResults_null_imp) <- saomResults_null_imp$effects$"effectName"
finalResults_null_imp$Estimate <- saomResults_null_imp$theta
finalResults_null_imp$SE <- saomResults_null_imp$se

#Odds Ratios
finalResults_null_imp$OddsRatio <- exp(finalResults_null_imp$Estimate)

#rounding
finalResults_null_imp_rounded <- as.data.frame(round(finalResults_null_imp, 3))
View(finalResults_null_imp_rounded)

#get a LateX ready table
stargazer(finalResults_null_imp_rounded, summary = FALSE, title="SAOM analysis, null imputation",
          rownames = TRUE)

#significance levels
finalResults_null_imp_rounded$sig <- NA
finalResults_null_imp_rounded$sig[(abs(finalResults_null_imp_rounded$Estimate) - abs(finalResults_null_imp_rounded$SE)*1.96)>0]<- "*"
finalResults_null_imp_rounded$sig[(abs(finalResults_null_imp_rounded$Estimate) - abs(finalResults_null_imp_rounded$SE)*2.58)>0]<- "**"
finalResults_null_imp_rounded$sig[(abs(finalResults_null_imp_rounded$Estimate) - abs(finalResults_null_imp_rounded$SE)*3.29)>0]<- "***"


