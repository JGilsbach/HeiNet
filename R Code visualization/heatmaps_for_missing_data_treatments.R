#heatmaps for missing data treatmenr/ imputation plots

#author: Judith Gilsbach

#packages needed
library(plotrix)
library(sysfonts)

#initialize fonts
font.add.google("Source Sans Pro")
windowsFonts(A = windowsFont("Source Sans Pro"),
             B = windowsFont("serif"))

#get Stataionary SAOM imputed data
load("StatSAOM_Imputation.RData")

#sum up matrices
ma_w3 <- Reduce('+', wave3imp_StatSAOM)
ma_w3_100 <- ma_w3
df_w3 <- as.data.frame(ma_w3_100)

#heatmap for wave 3: define colors
cellcol<-ma_w3_100
cellcol[cellcol== 50] <- "#597490" #observed values in blue
#works only because no value was imputed 50 times, otherwise change observed values to 100 beforehand
cellcol[cellcol>= 40 & cellcol <= 49] <-"#D0661C"
cellcol[cellcol>= 30 & cellcol <= 39] <-"#D98644"
cellcol[cellcol>= 20 & cellcol <= 29] <-"#E3A770"
cellcol[cellcol>= 10 & cellcol <= 19] <-"#EDC59F"
cellcol[cellcol>= 1 & cellcol <= 9] <-"#F6E3CF"
cellcol[cellcol== 0] <- "#FFFFFF"
  
testcol<-c("#597490","#D0661C","#D98644","#E3A770","#EDC59F","#F6E3CF","#FFFFFF")
col.labels<-c("Observed","Imputed min. 40 times","Imputed 30-39 times","Imputed 20-29 times","Imputed 10-19 times",
              "Imputed 1-9 times", "not observed or imputed") 

#missing actors wave 3  
x <- c(1,15,20,24,33,34,38,39,40)

#save plot
png(filename = "W3_StatSAOM.png",
    width = 8, height = 6, units = "in", pointsize = 12,
    bg = "transparent", res = 300)

par(mar=c(4,4,4,15))
par(family ="A")
color2D.matplot(df_w3, 
                show.values = F,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                na.color = "#FFFFFF",
                cellcolors=cellcol,
                show.legend=FALSE)

axis(3, at = seq_len(ncol(df_w3)) - 0.6,
     labels = names(df_w3), tick = FALSE, cex.axis = 0.6, mgp=c(0.1,0.1,0),las=2)
axis(2, at = seq_len(nrow(df_w3)) -0.5,
     labels = rev(rownames(df_w3)), col.axis = "black", tick = FALSE, las = 1, cex.axis = 0.6, mgp=c(0.1,0.1,0))
axis(2, at = c(42-0.5, 42-14.5,42-19.5,42-23.5,42-32.5,42-33.5,42-37.5,42-38.5,42-39.5),
     labels = x, col.axis = "#D0661C", tick = FALSE, las = 1, cex.axis = 0.6, mgp=c(0.1,0.1,0))

color.legend(45,15,50,30,rev(col.labels), rev(testcol), align = "rb",gradient="y")                        
                        
dev.off()  

#prepare data for plotting: Deletion of missing actors and null imputation
deletion_df <- as.data.frame(Bekannte_full_w3_ma)
deletion_df$`1` <- 2
deletion_df$`15` <- 2
deletion_df$`20` <- 2
deletion_df$`24` <- 2
deletion_df$`33` <- 2
deletion_df$`34` <- 2
deletion_df$`38` <- 2
deletion_df$`39` <- 2
deletion_df$`40` <- 2
deletion_ma <- as.matrix(deletion_df)
deletion_ma[is.na(deletion_ma)] <- 2


######## Plot with two panels: How null imputation and deletion method change the adjecency matrix #########
#color matrix for null imputation NA == 0
cellcol<- deletion_ma
cellcol[cellcol== 2] <- "#FFFFFF" #works only because no value was imputed 50 times
cellcol[cellcol== 1] <- "#597490"
cellcol[cellcol== 0] <- "#FFFFFF"

#cell color matrix for deletion of missing actors
cellcol2<- deletion_ma
cellcol2[cellcol2== 2] <- "#B5C1D3" #works only because no value was imputed 50 times
cellcol2[cellcol2== 1] <- "#597490"
cellcol2[cellcol2== 0] <- "#FFFFFF"
      
png(filename = "W3_Nullimp_Deletion_hochformat.png",
    width = 6, height = 12, units = "in", pointsize = 12,
    bg = "transparent", res = 300)

par(mar=c(2,2,2,2), mfrow=c(2,1))

par(family ="A")
color2D.matplot(deletion_ma, 
                show.values = F,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                na.color = "#FFFFFF",
                cellcolors=cellcol,
                show.legend=FALSE)

axis(3, at = seq_len(ncol(df_w3)) - 0.6,
     labels = names(df_w3), tick = FALSE, cex.axis = 0.6, mgp=c(0.1,0.1,0),las=2)
axis(2, at = seq_len(nrow(df_w3)) -0.5,
     labels = rev(rownames(df_w3)), col.axis = "black", tick = FALSE, las = 1, cex.axis = 0.6, mgp=c(0.1,0.1,0)) #rev(rownames(df_w3))
axis(2, at = c(42-0.5, 42-14.5,42-19.5,42-23.5,42-32.5,42-33.5,42-37.5,42-38.5,42-39.5),
     labels = x, col.axis = "#D0661C", tick = FALSE, las = 1, cex.axis = 0.6, mgp=c(0.1,0.1,0)) #rev(rownames(df_w3))



color2D.matplot(deletion_ma, 
                show.values = F,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                na.color = "#FFFFFF",
                cellcolors=cellcol2,
                show.legend=FALSE)

axis(3, at = seq_len(ncol(df_w3)) - 0.6,
     labels = names(df_w3), tick = FALSE, cex.axis = 0.6, mgp=c(0.1,0.1,0),las=2)
axis(2, at = seq_len(nrow(df_w3)) -0.5,
     labels = rev(rownames(df_w3)), col.axis = "black", tick = FALSE, las = 1, cex.axis = 0.6, mgp=c(0.1,0.1,0)) #rev(rownames(df_w3))
axis(2, at = c(42-0.5, 42-14.5,42-19.5,42-23.5,42-32.5,42-33.5,42-37.5,42-38.5,42-39.5),
     labels = x, col.axis = "#D0661C", tick = FALSE, las = 1, cex.axis = 0.6, mgp=c(0.1,0.1,0)) #rev(rownames(df_w3))

dev.off()  
                        