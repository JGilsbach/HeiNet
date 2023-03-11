##################### Plotting the SAOM results instead of tables #####################

#Author: Judith Gilsbach

library(tidyverse)
library(ggplot2)
library(RSiena)
library(knitr) #for rounding 

windowsFonts(A = windowsFont("Source Sans Pro"), #Gesis corporate font
             B = windowsFont("serif"))

#to avoid overload only the Stationary SAOM, null imputation and deletion method results are plotted

###load data: Stationary SAOM
load("StatSAOM_Imputation_Analysis_Results.RData")

###Deletion method
load("Deletion_SAOM_Analysis_Results.RData")

finalResults_deletion <- as.data.frame(matrix(,8,2))
names(finalResults_deletion) <- c("Estimate", "SE")
rownames(finalResults_deletion) <- saomResults_delete$effects$"effectName"
finalResults_deletion$Estimate <- saomResults_delete$theta
finalResults_deletion$SE <- saomResults_delete$se
finalResults_deletion_rounded <- as.data.frame(round(finalResults_deletion, 3))

#Odds Ratios
finalResults_deletion$OddsRatio <- exp(finalResults_deletion$Estimate)
finalResults_deletion_rounded <- as.data.frame(round(finalResults_deletion, 3))
View(finalResults_deletion_rounded)

### Null imputation
setwd("C:/Users/Judith Gilsbach/Desktop/Master/Modul 6 Masterarbeit/Saved Data/Test mit Null Imputation")
load("CompleteNullImp_SAOM_Analysis_Results.RData")

finalResults_null_imp <- as.data.frame(matrix(,8,2))
names(finalResults_null_imp) <- c("Estimate", "SE")
rownames(finalResults_null_imp) <- saomResults_null_imp$effects$"effectName"
finalResults_null_imp$Estimate <- saomResults_null_imp$theta
finalResults_null_imp$SE <- saomResults_null_imp$se


finalResults_null_imp_rounded <- as.data.frame(round(finalResults_null_imp, 3))

#Odds Ratios
finalResults_null_imp$OddsRatio <- exp(finalResults_null_imp$Estimate)
finalResults_null_imp_rounded <- as.data.frame(round(finalResults_null_imp, 3))
View(finalResults_null_imp_rounded)


####### Wrangling data for plotting ##########

#Data for SAOM results plot
save(finalResults_null_imp_rounded, finalResults_deletion_rounded,
     finalResults_SAOM_Int_rounded, finalResultsStatSAOM_rounded,
     file = "SAOM_Results_for_plotting.RData")

#combined data table
Est_null <- finalResults_null_imp_rounded$Estimate
Est_null <- Est_null[-7]
SE_null <- finalResults_null_imp_rounded$SE
SE_null <- SE_null[-7]
Est_deletion <- finalResults_deletion_rounded$Estimate
Est_deletion <- Est_deletion[-7]
SE_deletion <- finalResults_deletion_rounded$SE
SE_deletion <- SE_deletion[-7]
Est_stationary <- finalResultsStatSAOM_rounded$combinedEstimate
SE_stationary <- finalResultsStatSAOM_rounded$combinedSE

#create combined data.frame
SAOM_plotting_data <- data.frame(Est_null,
                                 SE_null,
                                 Est_deletion,
                                 SE_deletion,
                                 Est_stationary,
                                 SE_stationary)

row.names(SAOM_plotting_data) <- row.names(finalResultsStatSAOM_rounded)

SAOM_plotting_data$Parameter <- row.names(SAOM_plotting_data)

#save data
save(SAOM_plotting_data,
     file = "SAOM_Results_for_plotting_combined.RData")

#make separate datasets for rate constants and parameters
SAOM_plotting_rate_constants <- SAOM_plotting_data[1:2,]
SAOM_plotting_parameters <- SAOM_plotting_data[3:7,]
SAOM_plotting_parameters$Parameter <- c("density", "reciprocity", "GWESP", "outdegree activity","same tutorial")

#wrangling for parameters
SAOM_plotting_parameters_long <-  SAOM_plotting_parameters %>%
  pivot_longer(cols = !Parameter,
               names_to = "Imputation",
               values_to = "value") %>%
  separate(Imputation, c("value_kind", "method"), sep = "_") %>%
  pivot_wider(names_from = value_kind, values_from = value) %>%
  filter(!Parameter== "outdegree activity") %>% #remove results for outdegree activity
  mutate(num = 1:12) %>%
  mutate(method = recode(method, null = "null imputation", deletion = "deletion method", stationary = "stationary SAOM"))

#wrangling for rate constants
SAOM_plotting_rate_constants_long <-  SAOM_plotting_rate_constants %>%
  pivot_longer(cols = !Parameter,
               names_to = "Imputation",
               values_to = "value") %>%
  separate(Imputation, c("value_kind", "method"), sep = "_") %>%
  pivot_wider(names_from = value_kind, values_from = value) %>%
  select(-SE)

SAOM_plotting_rate_constants_long$Period<- NA
SAOM_plotting_rate_constants_long$Period[SAOM_plotting_rate_constants_long$Parameter =="constant Bekannte_full rate (period 1)"] <- "Period 1"
SAOM_plotting_rate_constants_long$Period[SAOM_plotting_rate_constants_long$Parameter =="constant Bekannte_full rate (period 2)"] <- "Period 2"

#Calculate 95% intervall for the rate constants of the multiply imputed values 
rate_constats_StatSAOM <- as.data.frame(matrix(,50,2))
names(rate_constats_StatSAOM) <- c("Period 1", "Period 2")
MIResults_StatSAOM[2,1]

z <- 1

for (x in 1:50){
  rate_constats_StatSAOM[x,1] <- MIResults_StatSAOM[1,z]
  rate_constats_StatSAOM[x,2] <- MIResults_StatSAOM[2,z]
  z <- z+2
}

#define confidence intervals 95%
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

confidence_interval(rate_constats_StatSAOM$`Period 1`,0.95)[[1]]
confidence_interval(rate_constats_StatSAOM$`Period 2`,0.95)
SAOM_plotting_rate_constants_long$lower <- NA
SAOM_plotting_rate_constants_long$upper <- NA
SAOM_plotting_rate_constants_long$lower[SAOM_plotting_rate_constants_long$method =="stationary" & 
                                       SAOM_plotting_rate_constants_long$Period == "Period 1"] <- confidence_interval(rate_constats_StatSAOM$`Period 1`,0.95)[[1]]
SAOM_plotting_rate_constants_long$lower[SAOM_plotting_rate_constants_long$method =="stationary" & 
                                       SAOM_plotting_rate_constants_long$Period == "Period 2"] <- confidence_interval(rate_constats_StatSAOM$`Period 2`,0.95)[[1]]
SAOM_plotting_rate_constants_long$upper[SAOM_plotting_rate_constants_long$method =="stationary" & 
                                          SAOM_plotting_rate_constants_long$Period == "Period 1"] <- confidence_interval(rate_constats_StatSAOM$`Period 1`,0.95)[[2]]
SAOM_plotting_rate_constants_long$upper[SAOM_plotting_rate_constants_long$method =="stationary" & 
                                          SAOM_plotting_rate_constants_long$Period == "Period 2"] <- confidence_interval(rate_constats_StatSAOM$`Period 2`,0.95)[[2]]


############################ Plots #####################

#only parameters
#wide x scale
setwd("C:/Users/Judith Gilsbach/Desktop/HeiNet")
png(filename = "Parameter Plot with legend above B.png",
    width = 20, height = 7, units = "in", pointsize = 12,
    bg = "transparent", res = 300)

ggplot(data=SAOM_plotting_parameters_long,
       aes(x = interaction(Parameter, num, lex.order = TRUE), y = Est, color = rev(method)), groups = 1)+ 
  geom_point(size = 5, shape = 19, position = position_dodge(w = 0.1)) +
  #annotate(geom = "text", x = seq_len(nrow(SAOM_plotting_parameters_long)), y = 3, label = SAOM_plotting_parameters_long$num, size = 0) +
  annotate(geom = "text", x = 2+3 * (0:3), y = -3.3, label = unique(SAOM_plotting_parameters_long$Parameter), size = 10) +
  coord_cartesian(ylim = c(-3, 3), expand = FALSE, clip = "off") +
  guides(color = guide_legend(override.aes = list(size=10, shape = 16),
                              title = "Imputation method"))+
  geom_errorbar(data=SAOM_plotting_parameters_long,aes(ymin=Est-SE,ymax=Est+SE),
                position = position_dodge(w = 0.1),
                width = 0.4, size = 1,
                show.legend = FALSE)+
  ylab("Estimate and SE") +
  xlab("Parameter")+
  scale_color_manual(values = c("#597490","#D0661C","#E3A770"))+
  geom_vline(xintercept = 3.5) +
  geom_vline(xintercept = 6.5) +
  geom_vline(xintercept = 9.5) +
  scale_x_discrete(breaks = seq(2, 11, by = 2)) +
  theme_bw(base_size = 30, base_family = "A") +
  theme(plot.margin = unit(c(1, 1, 3, 1), "lines"),
        axis.title.x = element_text(vjust=-4.5),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'top',
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"))


dev.off()    

#only rate constants
setwd("C:/Users/Judith Gilsbach/Desktop/HeiNet")
png(filename = "Rate Constant Plot with error bars.png",
    width = 5.5, height = 7, units = "in", pointsize = 12,
    bg = "transparent", res = 300)

ggplot(data=SAOM_plotting_rate_constants_long,aes(x = Period, y = Est, color = method))+ 
  geom_point(size = 5, shape = 19, position = position_dodge(w = 0.1)) +
  guides(color = guide_legend(override.aes = list(size=10, shape = 16),
                              title = "Imputation method"))+
  geom_errorbar(data=SAOM_plotting_rate_constants_long,aes(ymin=lower,ymax=upper), #add error bars for the multiply imputed model: Stationary SAOM
                position = position_dodge(w = 0.1),
                width = 1, size = 1,
                show.legend = FALSE)+
  ylab("Rate Constant \n (pooled if mdoel based)") +
  scale_color_manual(values = c("#D0661C","#E3A770","#597490","#B5C1D3"))+
  theme_minimal(base_size = 30, base_family = "A")+
  theme(legend.position = "none")

dev.off()

