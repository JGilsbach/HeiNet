#analysis with SAOM and with STERGM of HeiNet under deletion of missing actors

#for comparison with the multiply imputed analyses
#all actors who did not fill in all three surveys were deleted

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

#missings wave 1: 1,20,24 and 34
#missings wave 2: 1,21 and 40
#missings wave 3: 1, 15, 20, 24, 33, 34, 38, 39 and 40

#Delete missing vertices wave 1
Bekannte_full_w1_ma_delete <- Bekannte_full_w1_ma
Bekannte_full_w1_ma_delete <- Bekannte_full_w1_ma_delete[-c(1, 15, 20, 21, 24, 33, 34, 38, 39, 40),
                                                         -c(1, 15, 20, 21, 24, 33, 34, 38, 39, 40)]

#Delete missing vertices wave 2
Bekannte_full_w2_ma_delete <- Bekannte_full_w2_ma
Bekannte_full_w2_ma_delete <- Bekannte_full_w2_ma_delete[-c(1, 15, 20, 21, 24, 33, 34, 38, 39, 40),
                                                         -c(1, 15, 20, 21, 24, 33, 34, 38, 39, 40)]

#Delete missing vertices wave 3
Bekannte_full_w3_ma_delete <- Bekannte_full_w3_ma
Bekannte_full_w3_ma_delete <- Bekannte_full_w3_ma_delete[-c(1, 15, 20, 21, 24, 33, 34, 38, 39, 40),
                                                         -c(1, 15, 20, 21, 24, 33, 34, 38, 39, 40)]

#convert to network
Bekannte_full_net_w1_delete <- as.network(Bekannte_full_w1_ma_delete)
Bekannte_full_net_w2_delete <- as.network(Bekannte_full_w2_ma_delete)
Bekannte_full_net_w3_delete <- as.network(Bekannte_full_w3_ma_delete)

#delete attribute data for missing vertices
Attributes_K1_Deletion <- Attributes_K1[-c(1, 15, 20, 21, 24, 33, 34, 38, 39, 40),]

#Add Statistik attribute to all networks
set.vertex.attribute(Bekannte_full_net_w1_delete, 'Statistik', Attributes_K1_Deletion$Statistik_Zero)
set.vertex.attribute(Bekannte_full_net_w2_delete, 'Statistik', Attributes_K1_Deletion$Statistik_Zero)
set.vertex.attribute(Bekannte_full_net_w3_delete, 'Statistik', Attributes_K1_Deletion$Statistik_Zero)

############ STERGM analysis for the deletion method ################

#estimate STERGM_ Period 1: wave 1 and 2
Result_tergm_Deletion_w12 <- tergm(list(Bekannte_full_net_w1_delete, Bekannte_full_net_w2_delete)~
                                     Form(~ edges
                                          + mutual
                                          + gwesp(decay = log(2), fixed = TRUE)
                                          + odegree1.5() 
                                          + nodematch('Statistik'))+
                                     Diss(~ edges
                                          + mutual
                                          + nodematch('Statistik')),
                                   estimate="CMLE")

summary(Result_tergm_Deletion_w12)

#estimate STERGM_ Period 2: wave 2 and 3
Result_tergm_Deletion_w23 <- tergm(list(Bekannte_full_net_w2_delete, Bekannte_full_net_w3_delete)~
                                     Form(~ edges
                                          + mutual
                                          + gwesp(decay = log(2), fixed = TRUE)
                                          + odegree1.5() 
                                          + nodematch('Statistik'))+
                                     Diss(~ edges
                                          + mutual
                                          + nodematch('Statistik')),
                                   estimate="CMLE")

summary(Result_tergm_Deletion_w23)

#save results
save(Result_tergm_Deletion_w12, Result_tergm_Deletion_w23,
     Bekannte_full_net_w1_delete, Bekannte_full_net_w2_delete, Bekannte_full_net_w3_delete,
     file = "Deletion_TERGM_Analysis_Results.RData")


#save results in table: Period 1: wave 1 and 2
finalResults_TERGM_Deletion_w12 <- as.data.frame(matrix(,8,2))
names(finalResults_TERGM_Deletion_w12) <- c("Estimate", "SE")
rownames(finalResults_TERGM_Deletion_w12) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                               "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_Deletion_w12$Estimate <- summary(Result_tergm_Deletion_w12)$coefficients[1:8]
finalResults_TERGM_Deletion_w12$SE <- summary(Result_tergm_Deletion_w12)$asyse

#Odds Ratios
finalResults_TERGM_Deletion_w12$OddsRatio <- exp(finalResults_TERGM_Deletion_w12$Estimate)

#rounding
finalResults_TERGM_Deletion_rounded_w12 <- as.data.frame(round(finalResults_TERGM_Deletion_w12, 3))
View(finalResults_TERGM_Deletion_rounded_w12)

#save results in table: Period 2: wave 2 and 3
finalResults_TERGM_Deletion_w23 <- as.data.frame(matrix(,8,2))
names(finalResults_TERGM_Deletion_w23) <- c("Estimate", "SE")
rownames(finalResults_TERGM_Deletion_w23) <- c("Form~edges", "Form~mutual", "Form~gwesp.fixed.0.693147180559945", "Form~odegree1.5", "Form~nodematch.Statistik",
                                               "Diss~edges", "Diss~mutual", "Diss~nodematch.Statistik")
finalResults_TERGM_Deletion_w23$Estimate <- summary(Result_tergm_Deletion_w23)$coefficients[1:8]
finalResults_TERGM_Deletion_w23$SE <- summary(Result_tergm_Deletion_w23)$asyse


#Odds Ratios
finalResults_TERGM_Deletion_w23$OddsRatio <- exp(finalResults_TERGM_Deletion_w23$Estimate)

#rounding
finalResults_TERGM_Deletion_rounded_w23 <- as.data.frame(round(finalResults_TERGM_Deletion_w23, 3))
View(finalResults_TERGM_Deletion_rounded_w23)

#save both periods' results in one table
finalResults_TERGM_Deletion_rounded_w12_w23 <- finalResults_TERGM_Deletion_rounded_w12
finalResults_TERGM_Deletion_rounded_w12_w23$Estimate23 <- finalResults_TERGM_Deletion_rounded_w23$Estimate
finalResults_TERGM_Deletion_rounded_w12_w23$SE23 <- finalResults_TERGM_Deletion_rounded_w23$SE
finalResults_TERGM_Deletion_rounded_w12_w23$OddsRatio23 <- finalResults_TERGM_Deletion_rounded_w23$OddsRatio
View(finalResults_TERGM_Deletion_rounded_w12_w23)
finalResults_TERGM_Deletion_rounded_w12_w23 <- as.data.frame(finalResults_TERGM_Deletion_rounded_w12_w23)

#get LaTeX ready table
stargazer(finalResults_TERGM_Deletion_rounded_w12_w23, summary = FALSE, title="TERGM analysis, missing actors deleted",
          rownames = TRUE, column.labels= c("Period 1", "Period 2"), column.separate = c(3,3))

#significance levels
finalResults_TERGM_Deletion_rounded_w12_w23$sig_w12 <- NA
finalResults_TERGM_Deletion_rounded_w12_w23$sig_w12[(abs(finalResults_TERGM_Deletion_rounded_w12_w23$Estimate) - abs(finalResults_TERGM_Deletion_rounded_w12_w23$SE)*1.96)>0]<- "*"
finalResults_TERGM_Deletion_rounded_w12_w23$sig_w12[(abs(finalResults_TERGM_Deletion_rounded_w12_w23$Estimate) - abs(finalResults_TERGM_Deletion_rounded_w12_w23$SE)*2.58)>0]<- "**"
finalResults_TERGM_Deletion_rounded_w12_w23$sig_w12[(abs(finalResults_TERGM_Deletion_rounded_w12_w23$Estimate) - abs(finalResults_TERGM_Deletion_rounded_w12_w23$SE)*3.29)>0]<- "***"

finalResults_TERGM_Deletion_rounded_w12_w23$sig_w23 <- NA
finalResults_TERGM_Deletion_rounded_w12_w23$sig_w23[(abs(finalResults_TERGM_Deletion_rounded_w12_w23$Estimate23) - abs(finalResults_TERGM_Deletion_rounded_w12_w23$SE23)*1.96)>0]<- "*"
finalResults_TERGM_Deletion_rounded_w12_w23$sig_w23[(abs(finalResults_TERGM_Deletion_rounded_w12_w23$Estimate23) - abs(finalResults_TERGM_Deletion_rounded_w12_w23$SE23)*2.58)>0]<- "**"
finalResults_TERGM_Deletion_rounded_w12_w23$sig_w23[(abs(finalResults_TERGM_Deletion_rounded_w12_w23$Estimate23) - abs(finalResults_TERGM_Deletion_rounded_w12_w23$SE23)*3.29)>0]<- "***"


############ SAOM analysis for the deletion method ################

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

#prepare for SAOM
Statistik_v <- as.vector(Attributes_K1_Deletion$Statistik)
Statistik  <- coCovar(as.vector(Statistik_v),centered=FALSE)

Bekannte_full <- sienaDependent(array(c(Bekannte_full_w1_ma_delete, Bekannte_full_w2_ma_delete,
                                        Bekannte_full_w3_ma_delete), dim = c( 32, 32, 3)))

Data  <- sienaDataCreate(Bekannte_full, Statistik)
effectsData <- getEffects(Data)
effectsData <- includeEffects(effectsData,
                              gwespFF, recip, density, outAct, outIso) #outIso (inverse outTrunc for p = 1) = positive coefficient means that there is a positive tendency toward outdegree = 0
effectsData <- includeEffects(effectsData,  sameX,
                              interaction1 =  "Statistik")

estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 312,
                                           n3 = 1000, maxlike = FALSE, cond = FALSE,
                                           lessMem = FALSE)

saomResults_delete <- siena07ToConvergence(alg = estimation.options,
                                           dat = Data, eff = effectsData,
                                           threshold = 0.25)

saomResults_delete

#colliearity check
summary(saomResults_delete)

#check for autocorrelation
saomResults_delete$ac

#goodness of fit
saom.results.gof.result_Delete <- list()
saom.results.gof.result_Delete[[1]] <- plot(sienaGOF(saomResults_delete, varName="Bekannte_full", OutdegreeDistribution))
saom.results.gof.result_Delete[[2]] <- plot(sienaGOF(saomResults_delete, varName="Bekannte_full", IndegreeDistribution))
saom.results.gof.result_Delete[[1]] 
saom.results.gof.result_Delete[[2]] 

#save results in data.frame
finalResults_deletion <- as.data.frame(matrix(,8,2))
names(finalResults_deletion) <- c("Estimate", "SE")
rownames(finalResults_deletion) <- saomResults_delete$effects$"effectName"
finalResults_deletion$Estimate <- saomResults_delete$theta
finalResults_deletion$SE <- saomResults_delete$se

#Odds Ratio
finalResults_deletion$OddsRatio <- exp(finalResults_deletion$Estimate)

#rounding
finalResults_deletion_rounded <- as.data.frame(round(finalResults_deletion, 3))

#save results
save(saomResults_delete, finalResults_deletion_rounded,
     Bekannte_full_net_w1_delete, Bekannte_full_net_w2_delete, Bekannte_full_net_w3_delete, saom.results.gof.result_Delete, 
     file = "Deletion_SAOM_Analysis_Results.RData")

#get LaTeX ready table
stargazer(finalResults_deletion_rounded, summary = FALSE, title="SAOM analysis, deletion of missing actors",
          rownames = TRUE)

#significance levels
finalResults_deletion_rounded$sig <- NA
finalResults_deletion_rounded$sig[(abs(finalResults_deletion_rounded$Estimate) - abs(finalResults_deletion_rounded$SE)*1.96)>0]<- "*"
finalResults_deletion_rounded$sig[(abs(finalResults_deletion_rounded$Estimate) - abs(finalResults_deletion_rounded$SE)*2.58)>0]<- "**"
finalResults_deletion_rounded$sig[(abs(finalResults_deletion_rounded$Estimate) - abs(finalResults_deletion_rounded$SE)*3.29)>0]<- "***"
