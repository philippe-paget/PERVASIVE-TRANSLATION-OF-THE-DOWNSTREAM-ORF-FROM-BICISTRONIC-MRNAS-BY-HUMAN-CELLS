### CODOVIREVOL - E6E7 project

### Flow cytometry analyses: GATING STEPS

#Put all your csv files in to the same folder
#Check that they all have the same name template 

###########
## SETUP ##
###########

## Load packages
library(ape) 
library(ade4)
library(vegan)
library(cluster)
library(MASS)
library(ggplot2)
library(stringr)
library(data.table)
library(glm2)
library(spatialEco) 
library(tidyverse)
library(reshape2)
library(gplots)
library(plyr)

## Define input data
#setwd("_______") # Put here the location of the folder where all your files are

setwd("~/XXX/")

temp = list.files(pattern="*.csv")
YOURDATA = lapply(temp, read.csv)
names(YOURDATA) <-(temp)
names(YOURDATA)
length(YOURDATA) #should be equal to 60

#changer le nom des colonnes
canard<-function(i){
  colnames(i)<-c("FSC.H", "SSC.H", "GFP.H", "FSC.A", "SSC.A", "GFP.A", "Width", "Time", "SSC.HoverSSC.A")
  return(i)
}
YOURDATA<-lapply(YOURDATA,canard)

Count.events <- function(c){      # here we count the number of events in each file
  nrow(c)}
YOURDATA.events <- lapply(YOURDATA, Count.events)
YOURDATA.events # should be equal or close to 50000 (number of events) 

## Preliminary plots
#####################

for (i in 1:xxx) {          
  samplename = names(YOURDATA)[i]
  par(mfcol = c(1, 3), mar = c(7, 5, 4, 2), oma = c(0, 0, 3, 0), mex=1)
  sample_test <- get(samplename,YOURDATA.1)
  plot(density(sample_test$FSC.H),
       main="FSC-H density plot",
       xlab="FSC-H (size)",
       ylab="Density of events")
  plot(density(as.numeric(sample_test$SSC.H)),
       main="SSC-H density plot",
       xlab="SSC-H (granulosity / complexity)",
       ylab="Density of events")
  plot(sample_test$SSC.H,sample_test$FSC.H,
       main="Cell populations",
       xlab="SSC-H (granulosity / complexity)",
       ylab="FSC-H (size)",xlim=c(0,2500000), ylim=c(0,10000000), pch=3,cex = 0.1)
  mtext(paste("Raw data - ",samplename,sep=""), side = 3, outer = TRUE, cex = 1)
}

                ############
                ## GATING ##
                ############

### STEP 1: Remove outliers
###########################

## Define function
Nooutliers <- function (i){
  i1<-subset(i, i$FSC.H == (winsorize(i$FSC.H, p=c(0.005, 0.995))) & i$SSC.H == (winsorize(i$SSC.H, p=c(0.005, 0.995))) 
       & i$FSC.A == (winsorize(i$FSC.A, p=c(0.005, 0.995))) & i$SSC.A == (winsorize(i$SSC.A, p=c(0.005, 0.995))))
}


## Apply function on list
YOURDATA.W <- lapply(YOURDATA, Nooutliers)
length(YOURDATA.W) # Should still be equal to 185 #553

## Count remaining events
YOURDATA.W.events <- lapply(YOURDATA.W, Count.events)
YOURDATA.W.events
YOURDATA.events <- lapply(YOURDATA.1, Count.events)
YOURDATA.events
YOURDATA.W.ratio <- as.numeric(YOURDATA.W.events) / as.numeric(YOURDATA.events)
YOURDATA.W.ratio # calculate the proportion of remaining events after outlier filtering
max(YOURDATA.W.ratio) 
min(YOURDATA.W.ratio) 

## Plots as pdf output files  # still not working for me...

for (i in 1:xxx) {
  samplename = names(YOURDATA.W)[i]
  par(mfcol = c(1, 3), mar = c(7, 5, 4, 2), oma = c(0, 0, 3, 0), mex=1)
  sample_test <- get(samplename,YOURDATA.W)
  plot(density(sample_test$FSC.H),
       main="FSC-H density plot",
       xlab="FSC-H (size)",
       ylab="Density of events")
  plot(density(as.numeric(sample_test$SSC.H)),
       main="SSC-H density plot",
       xlab="SSC-H (granulosity / complexity)",
       ylab="Density of events")
  plot(sample_test$SSC.H,sample_test$FSC.H,
       main="Cell populations",
       xlab="SSC-H (granulosity / complexity)",
       ylab="FSC-H (size)",xlim=c(0,2500000), ylim=c(0,10000000), pch=3,cex = 0.1)
  mtext(paste("No outliers - ",samplename,sep=""), side = 3, outer = TRUE, cex = 1)
}

for (i in 1:xxx) {
  samplename = names(YOURDATA.W)[i]
  par(mfcol = c(1, 1), mar = c(7, 5, 4, 2), oma = c(0, 0, 3, 0), mex=1)
  sample_test <- get(samplename,YOURDATA.W)
  
  x <- densCols(sample_test$SSC.H,sample_test$FSC.H) 
  sample_test$dens <- col2rgb(x)[1,] + 1L
  cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(100) #(50) #low number, hot colors #high number, cold colors
  sample_test$col <- cols[sample_test$dens]
  plot(sample_test$FSC.H~sample_test$SSC.H, data=sample_test[order(sample_test$dens),], pch=20, col=col, cex=1,
       main="Cell populations - Without outliers",
       xlab="SSC-H (granulosity/complexity)",
       ylab="FSC-H (size)")
  
  mtext(paste("No outliers - ",samplename,sep=""), side = 3, outer = TRUE, cex = 2)
}

### STEP 2: Remove debris
#########################

## Define functions
Remove.Debris <- function(i){
  dSscH <- density(i$SSC.H)
  dFscH <- density(i$FSC.H)
  dSscHT <- as.data.frame(cbind(dSscH$x, dSscH$y))
  dFscHT <- as.data.frame(cbind(dFscH$x, dFscH$y)) 
  # make tables with abscissa in first column (V1) and ordinate in second column (V2)
  minFy <- (local.min.max(dFscH$y, plot = FALSE))$minima[1] 
  maxFy <- max((local.min.max(dFscH$y, plot = FALSE))$maxima)
  minSy <- (local.min.max(dSscH$y, plot = FALSE))$minima[1]
  maxSy  <- max((local.min.max(dSscH$y, plot = FALSE))$maxima)
  # determine the minimum and maximum densities of FSC-H & SSC-H (deeps & peaks of the curve)
  minFx <- subset(dFscHT, dFscHT$V2 == minFy)$V1
  maxFx <- (subset(dFscHT, dFscHT$V2 == maxFy))$V1 
  minSx <- subset(dSscHT, dSscHT$V2 == minSy)$V1
  maxSx <- (subset(dSscHT, dSscHT$V2 == maxSy))$V1
  # determine the abscissa corresponding to the minimum and maximum densities
  j <- i[-which(i$FSC.H < minFx),]
  h <- j[-which(j$SSC.H < minSx),]
}

#recuperer les valeurs minimales
minimaF <- function(x){
  dFscH<-density(x$FSC.H)
  dFscHT <- as.data.frame(cbind(dFscH$x, dFscH$y))
  minFy <- (local.min.max(dFscH$y, plot = FALSE))$minima[1]
  minF2 <- subset(dFscHT, dFscHT$V2 == minFy)
  minF3 <- minF2$V1
}

minimaS <- function(x){
  dSscH<-density(x$SSC.H)
  dSscHT <- as.data.frame(cbind(dSscH$x, dSscH$y))
  minSy <- (local.min.max(dSscH$y, plot = FALSE))$minima[1]
  minS2 <- subset(dSscHT, dSscHT$V2 == minSy)
  minS3 <- minS2$V1  
}

#Calculer les moyennes des vecteurs
DATAMINIMAF <- lapply(YOURDATA.W, minimaF)
DATAMINIMAS <- lapply(YOURDATA.W , minimaS)
papillon = mean(unlist(lapply(DATAMINIMAF,mean)),na.rm=T)
libellule = mean(unlist(lapply(DATAMINIMAS,mean)),na.rm=T)

#enlever les valeurs inferieures aux moyennes

## Apply function on list
YOURDATA.D0 <- lapply(YOURDATA.W, Remove.Debris)

## Plots

for (i in 1:xxx) {
  samplename = names(YOURDATA.D0)[i]
  par(mfcol = c(1, 3), mar = c(7, 5, 4, 2), oma = c(0, 0, 3, 0), mex=1)
  sample_test <- get(samplename,YOURDATA.D0)
  plot(density(sample_test$FSC.H),
       main="FSC-H density plot",
       xlab="FSC-H (size)",
       ylab="Density of events")
  plot(density(as.numeric(sample_test$SSC.H)),
       main="SSC-H density plot",
       xlab="SSC-H (granulosity / complexity)",
       ylab="Density of events")
  plot(sample_test$SSC.H,sample_test$FSC.H,
       main="Cell populations",
       xlab="SSC-H (granulosity/complexity)",
       ylab="FSC-H (size)")
  mtext(paste("No debris - ",samplename,sep=""), side = 3, outer = TRUE, cex = 1)
}


## Plots in color

for (i in 1:xxx) {
  samplename = names(YOURDATA.D0)[i]
  par(mfcol = c(1, 1), mar = c(7, 5, 4, 2), oma = c(0, 0, 3, 0), mex=1)
  sample_test <- get(samplename,YOURDATA.D0)
  
  x <- densCols(sample_test$SSC.H,sample_test$FSC.H) 
  sample_test$dens <- col2rgb(x)[1,] + 1L
  cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(100) #(50) #low number, hot colors #high number, cold colors
  sample_test$col <- cols[sample_test$dens]
  plot(sample_test$FSC.H~sample_test$SSC.H, data=sample_test[order(sample_test$dens),], pch=20, col=col, cex=1,
       main="Cell populations - Without debris",
       xlab="SSC-H (granulosity/complexity)",
       ylab="FSC-H (size)")
  
  mtext(paste("No debris - ",samplename,sep=""), side = 3, outer = TRUE, cex = 2)
}

## Estimate the number of removed events
Count.events <- function(c){
  nrow(c)}
YOURDATA.D0.events <- lapply(YOURDATA.D0, Count.events)
YOURDATA.D0.events
YOURDATA.W.events <- lapply(YOURDATA.W, Count.events)
YOURDATA.W.events
YOURDATA.D0.prop <- as.numeric(YOURDATA.D0.events)/as.numeric(YOURDATA.W.events)
YOURDATA.D0.prop
length (YOURDATA.D0.prop)
YOURDATA.D0.prop.table <- as.data.frame(cbind(names(YOURDATA.D0), YOURDATA.D0.prop))
names(YOURDATA.D0.prop.table) <- c("sample","conserved_prop")
head(YOURDATA.D0.prop.table)
YOURDATA.many.debris <- subset(YOURDATA.D0.prop.table, as.numeric(as.character(YOURDATA.D0.prop.table$conserved_prop))<=0.80)
YOURDATA.many.debris
# check if samples with more than 20% of filtered-out events display good graphical outputs (see visualisation of the data above)


### STEP 3: Remove doublets
###########################

## Define function 1 (standard)

remove.doublets <- function(D) {
  
  D <- subset(D,D$SSC.H>=0 & D$SSC.A>=0)
  D$`SSC.H/SSC.A` <- D$SSC.H/D$SSC.A
  D <- subset(D,D$`SSC.H/SSC.A` == (winsorize(D$`SSC.H/SSC.A`,p=c(0.005, 0.995))))
  # ratio of cell size/cell signal
  dW <- density(D$`SSC.H/SSC.A`)
  # distribution of this ratio in the sample: plot(dW) display a main peak and a smaller peak on its left
  dWT <- as.data.frame(cbind(dW$x, dW$y))
  # make a table with the coordinates V1 = abscissa and V2 = ordinates
  minmax <- c((local.min.max(dWT$V2, plot = FALSE))$minima, (local.min.max(dWT$V2, plot = FALSE))$maxima)
  # calculation of local min and max
  minmaxT <- dWT[dWT$V2 %in% minmax,]
  # subset the coordinate table with only coordinates corresponding to min or max local (V2 = ordinates)
  maxpoint <- max(minmaxT$V2)
  maxpointx<- (subset(minmaxT, minmaxT$V2 == maxpoint))$V1
  # the max of the local max and its abscissa
  maxpos <- match(maxpointx, minmaxT$V1)
  # position of the max abscissa in the vector of all abscissa (or which ligne in the first column)
  lentlim = minmaxT$V1[maxpos-1]
  # determine the local minimum at the left of the the highest density peak / conditions are necessary because the peak is double
  D2 <- subset(D, (D$`SSC.H/SSC.A` > lentlim))
  # window to select
  
}

## Apply function
YOURDATA.CLEAN.IND <- lapply(YOURDATA.D0, remove.doublets)
length(YOURDATA.CLEAN.IND)
YOURDATA.CLEAN.IND_1 <- subset(YOURDATA.CLEAN.IND,names(YOURDATA.CLEAN.IND)!="xxx")
length (YOURDATA.CLEAN.IND_1)

## Define function 2 (for the exceptions defined above)
YOURDATA.D0_2 <- subset(YOURDATA.D0,names(YOURDATA.D0)=="xxx") # put the name of the exceptions
length (YOURDATA.D0_2)

remove.doublets <- function(D) {
  
  D <- subset(D,D$SSC.H>=0 & D$SSC.A>=0)
  D$`SSC.H/SSC.A` <- D$SSC.H/D$SSC.A
  D <- subset(D,D$`SSC.H/SSC.A` == (winsorize(D$`SSC.H/SSC.A`,p=c(0.005, 0.995))))
  # ratio of cell size/cell signal
  dW <- density(D$`SSC.H/SSC.A`)
  # distribution of this ratio in the sample: plot(dW) display a main peak and a smaller peak on its left
  dWT <- as.data.frame(cbind(dW$x, dW$y))
  # make a table with the coordinates V1 = abscissa and V2 = ordinates
  minmax <- c((local.min.max(dWT$V2, plot = FALSE))$minima, (local.min.max(dWT$V2, plot = FALSE))$maxima)
  # calculation of local min and max
  minmaxT <- dWT[dWT$V2 %in% minmax,]
  # subset the coordinate table with only coordinates corresponding to min or max local (V2 = ordinates)
  maxpoint <- max(minmaxT$V2)
  maxpointx<- (subset(minmaxT, minmaxT$V2 == maxpoint))$V1
  # the max of the local max and its abscissa
  maxpos <- match(maxpointx, minmaxT$V1)
  # position of the max abscissa in the vector of all abscissa (or which ligne in the first column)

if ((maxpos-3) <= 0) {lentlim = minmaxT$V1[maxpos-1]} else {
  if (minmaxT$V1[maxpos-1] > minmaxT$V1[maxpos-3]) {lentlim = minmaxT$V1[maxpos-3]} else {
    lentlim = minmaxT$V1[maxpos-1]}
  }
  # determine the local minimum at the left of the the highest density peak / conditions are necessary because the peak is double
  D2 <- subset(D, (D$`SSC.H/SSC.A` > lentlim))
  # window to select

}

## Apply function
YOURDATA.CLEAN.IND_2 <- lapply(YOURDATA.D0_2, remove.doublets)
length(YOURDATA.CLEAN.IND_2)

## Merge
YOURDATA.CLEAN.IND = c(YOURDATA.CLEAN.IND_1, YOURDATA.CLEAN.IND_2)
names(YOURDATA.CLEAN.IND)
save(YOURDATA.CLEAN.IND, file = "R8_2022_E6E7")
## Visualize
# Names of exceptions: "xxx|"xxx"|"xxx"|"xxx")
sample_test=xxx
names(sample_test)
head(sample_test)
plot(sample_test$SSC.A,sample_test$SSC.H) #after
plot(density(sample_test$SSC.H/sample_test$SSC.A))

D=YOURDATA.D0_2$"xxx"
D <- subset(D,D$SSC.H>=0 & D$SSC.A>=0)
D$`SSC.H/SSC.A` <- D$SSC.H/D$SSC.A
D <- subset(D,D$`SSC.H/SSC.A` == (winsorize(D$`SSC.H/SSC.A`,p=c(0.005, 0.995))))
plot(D$SSC.A,D$SSC.H) # before
plot(density(D$SSC.H/D$SSC.A))

for (i in 1:xxx) {
  # general parameters
  par(mfcol = c(1, 2), mar = c(7, 5, 4, 2), oma = c(0, 0, 3, 0), mex=1)
  samplename0 = names(YOURDATA.D0)[i]
  sample_test <- get(samplename0,YOURDATA.CLEAN.IND)
  sample_test0 <- get(samplename0,YOURDATA.D0)
  #samplename = names(YOURDATA.CLEAN.IND_2)[i] # for the 5 exceptions
  #samplename0 = names(YOURDATA.D0_2)[i] # for the 5 exceptions
  #sample_test <- get(samplename,YOURDATA.CLEAN.IND_2) # for the 5 exceptions
  #sample_test0 <- get(samplename0,YOURDATA.D0_2) # for the 5 exceptions
  xmax <- max(sample_test0$SSC.A) + 1000
  xmin <- min(sample_test0$SSC.A) 
  ymax <- max(sample_test0$SSC.H) + 1000
  ymin <- min(sample_test0$SSC.H) 
  # black & white
  plot(sample_test0$SSC.A,sample_test0$SSC.H,
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),
       cex=0.5,
       main="Before doublet gating",
       xlab="SSC-A (event signal)",
       ylab="SSC-H (event size)")
  plot(sample_test$SSC.A,sample_test$SSC.H,
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),
       cex=0.5,
       main="After doublet gating",
       xlab="SSC-A (event signal)",
       ylab="SSC-H (event size)")
  mtext(paste("With/Without doublets comparison - ", samplename0, sep=""), side = 3, outer = TRUE, cex = 1)
}

# colors 

for (i in 1:xxx) {
  # general parameters
  par(mfcol = c(1, 2), mar = c(7, 5, 4, 2), oma = c(0, 0, 3, 0), mex=1)
  samplename0 = names(YOURDATA.D0)[i]
  sample_test <- get(samplename0,YOURDATA.CLEAN.IND)
  sample_test0 <- get(samplename0,YOURDATA.D0)
  #samplename = names(YOURDATA.CLEAN.IND_2)[i] # for the 5 exceptions
  #samplename0 = names(YOURDATA.D0_2)[i] # for the 5 exceptions
  #sample_test <- get(samplename,YOURDATA.CLEAN.IND_2) # for the 5 exceptions
  #sample_test0 <- get(samplename0,YOURDATA.D0_2) # for the 5 exceptions
  xmax <- max(sample_test0$SSC.A) + 1000
  xmin <- min(sample_test0$SSC.A) 
  ymax <- max(sample_test0$SSC.H) + 1000
  ymin <- min(sample_test0$SSC.H) 
 # color plot before
x <- densCols(sample_test0$SSC.A,sample_test0$SSC.H) 
sample_test0$dens <- col2rgb(x)[1,] + 1L
cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(100) #(50) #low number, hot colors #high number, cold colors
sample_test0$col <- cols[sample_test0$dens]
plot(sample_test0$SSC.H~sample_test0$SSC.A, data=sample_test0[order(sample_test0$dens),],
     pch=20, col=col, cex=0.5,
     xlim=c(xmin,xmax),ylim=c(ymin,ymax),
     main="Before doublet gating",
     xlab="SSC-A (event signal)",
     ylab="SSC-H (event size)")
  # color plot after
x <- densCols(sample_test$SSC.A,sample_test$SSC.H) 
sample_test$dens <- col2rgb(x)[1,] + 1L
cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(100) #(50) #low number, hot colors #high number, cold colors
sample_test$col <- cols[sample_test$dens]
plot(sample_test$SSC.H~sample_test$SSC.A, data=sample_test[order(sample_test$dens),],
     pch=20, col=col, cex=0.5,
     xlim=c(xmin,xmax),ylim=c(ymin,ymax),
     main="After doublet gating",
     xlab="SSC-A (event signal)",
     ylab="SSC-H (event size)")
mtext(paste("With/Without doublets comparison - ", samplename0, sep=""), side = 3, outer = TRUE, cex = 1)

}


# density curves 
for (i in 1:xxx) {
  par(mfcol = c(1, 2), mar = c(7, 5, 4, 2), oma = c(0, 0, 3, 0), mex=1)
  samplename0 = names(YOURDATA.D0)[i]
  sample_test <- get(samplename0,YOURDATA.CLEAN.IND)
  sample_test0 <- get(samplename0,YOURDATA.D0)
  #samplename = names(YOURDATA.CLEAN.IND_2)[i] # for the 5 exceptions
  #samplename0 = names(YOURDATA.D0_2)[i] # for the 5 exceptions
  #sample_test <- get(samplename,YOURDATA.CLEAN.IND_2) # for the 5 exceptions
  #sample_test0 <- get(samplename0,YOURDATA.D0_2) # for the 5 exceptions
  sample_test0 <- subset(sample_test0,sample_test0$SSC.H>=0 & sample_test0$SSC.A>=0)
  sample_test0$`SSC.H/SSC.A` <- sample_test0$SSC.H/sample_test0$SSC.A
  sample_test0 <- subset(sample_test0,sample_test0$`SSC.H/SSC.A` == (winsorize(sample_test0$`SSC.H/SSC.A`,p=c(0.005, 0.995))))
  dW0 <- density(sample_test0$`SSC.H/SSC.A`)
  plot(dW0,xlim=c(1,3.5),cex=0.2,main="Before gating")
  sample_test <- subset(sample_test,sample_test$SSC.H>=0 & sample_test$SSC.A>=0)
  sample_test$`SSC.H/SSC.A` <- sample_test$SSC.H/sample_test$SSC.A
  sample_test <- subset(sample_test,sample_test$`SSC.H/SSC.A` == (winsorize(sample_test$`SSC.H/SSC.A`,p=c(0.005, 0.995))))
  dW <- density(sample_test$`SSC.H/SSC.A`)
  plot(dW,xlim=c(1,3.5),cex=0.2,main="After gating")
  mtext(paste("With/Without doublets comparison - ", samplename0, sep=""), side = 3, outer = TRUE, cex = 1)
}  
  
  ### SUMMARY : Estimate the number of removed events 

## TABLE COUNT
# Gating steps (for reminder):
  #0. Raw data: YOURDATA | YOURDATA.1
  #1. Remove outliers: YOURDATA.W
  #2. Remove debris: YOURDATA.D0
  #3. Remove doublets: YOURDATA.CLEAN.IND 

IncrementalTable<-c()
for (i in 1:xxx) {
  samplename0 = names(YOURDATA.D0)[i]
  sample_doublets <- get(samplename0,YOURDATA.CLEAN.IND)
  sample_debris <- get(samplename0,YOURDATA.D0)
  sample_outliers <- get(samplename0,YOURDATA.W)
  sample_raw <- get(samplename0,YOURDATA)
  
  sample_doublets_num <- nrow(sample_doublets)
  sample_debris_num <- nrow(sample_debris)
  sample_outliers_num <- nrow(sample_outliers)
  sample_raw_num <- nrow(sample_raw)
  
  proportion_of_outliers <- round((1-as.numeric(sample_outliers_num)/as.numeric(sample_raw_num))*100,digits=2)
  proportion_of_debris <- round((1-as.numeric(sample_debris_num)/as.numeric(sample_outliers_num))*100,digits=2)
  proportion_of_doublets <- round((1-(as.numeric(sample_doublets_num)/as.numeric(sample_debris_num)))*100,digits=2)
  proportion_of_conserved_events <- round((as.numeric(sample_doublets_num)/as.numeric(sample_raw_num))*100,digits=2)

  IncrementalTable <- rbind(IncrementalTable,c(samplename0,sample_raw_num,sample_outliers_num,sample_debris_num,sample_doublets_num,proportion_of_outliers,proportion_of_debris,proportion_of_doublets,proportion_of_conserved_events))
}

rownames(IncrementalTable) <- IncrementalTable[,1]
colnames(IncrementalTable) <- c("sample","raw_events","no_outliers","no_debris","no_doublets","proportion_of_outliers","proportion_of_debris","proportion_of_doublets","proportion_of_conserved_events")
head(IncrementalTable)
IncrementalTable

mypath="~/xxx/"
write.table(IncrementalTable,file=paste(mypath,"table_count",".csv",sep=""),col.names = TRUE,row.names = FALSE)

## PICS PER SAMPLE

mypath="~/xxx/"
pdf(file=paste(mypath,"gating_summary_per_sample",".pdf",sep=""),width=8.27,height = 11.69)

for (i in 1:xxx) {
  par(mfrow = c(4, 3), mar = c(7, 5, 4, 2), oma = c(0, 0, 3, 0), mex=1)

# raw data curves & graph
  samplename0 = names(YOURDATA)[i]
    sample_test <- get(samplename0,YOURDATA)
    plot(density(sample_test$FSC.H),
         main="Size distribution - raw data",
         xlab="FSC-H (size)",
         ylab="Density of events")
    plot(density(as.numeric(sample_test$SSC.H)),
         main="Complexity distribution - raw data",
         xlab="SSC-H (complexity)",
         ylab="Density of events")
    x <- densCols(sample_test$SSC.H,sample_test$FSC.H) 
    sample_test$dens <- col2rgb(x)[1,] + 1L
    cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(100)
    sample_test$col <- cols[sample_test$dens]
    plot(sample_test$FSC.H~sample_test$SSC.H, data=sample_test[order(sample_test$dens),], pch=20, col=col, cex=0.5,
         main="Cell populations - raw data",
         xlab="SSC-H (complexity)",
         ylab="FSC-H (size)")

# outliers removal curves & graph
    samplename0 = names(YOURDATA)[i]
    sample_test <- get(samplename0,YOURDATA.W)
    plot(density(sample_test$FSC.H),
         main="Size distribution - no outliers",
         xlab="FSC-H (size)",
         ylab="Density of events")
    plot(density(as.numeric(sample_test$SSC.H)),
         main="Complexity distribution - no outliers",
         xlab="SSC-H (complexity)",
         ylab="Density of events")
    x <- densCols(sample_test$SSC.H,sample_test$FSC.H) 
    sample_test$dens <- col2rgb(x)[1,] + 1L
    cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(100)
    sample_test$col <- cols[sample_test$dens]
    plot(sample_test$FSC.H~sample_test$SSC.H, data=sample_test[order(sample_test$dens),], pch=20, col=col, cex=0.5,
         main="Cell populations - no outliers",
         xlab="SSC-H (complexity)",
         ylab="FSC-H (size)")
    
# debris removal curves & graph
    samplename0 = names(YOURDATA)[i]
    sample_test <- get(samplename0,YOURDATA.D0)
    plot(density(sample_test$FSC.H),
         main="Size distribution - no debris",
         xlab="FSC-H (size)",
         ylab="Density of events")
    plot(density(as.numeric(sample_test$SSC.H)),
         main="Complexity distribution - no debris",
         xlab="SSC-H (complexity)",
         ylab="Density of events")
    x <- densCols(sample_test$SSC.H,sample_test$FSC.H) 
    sample_test$dens <- col2rgb(x)[1,] + 1L
    cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(100)
    sample_test$col <- cols[sample_test$dens]
    plot(sample_test$FSC.H~sample_test$SSC.H, data=sample_test[order(sample_test$dens),], pch=20, col=col, cex=0.5,
         main="Cell populations - no debris",
         xlab="SSC-H (complexity)",
         ylab="FSC-H (size)")
    
# doublet removal curves
  samplename0 = names(YOURDATA)[i]
  sample_test <- get(samplename0,YOURDATA.CLEAN.IND)
  sample_test0 <- get(samplename0,YOURDATA.D0)
  # before
  sample_test1 <- subset(sample_test0,sample_test0$SSC.H>=0 & sample_test0$SSC.A>=0)
  sample_test1$`SSC.H/SSC.A` <- sample_test1$SSC.H/sample_test1$SSC.A
  sample_test1 <- subset(sample_test1,sample_test1$`SSC.H/SSC.A` == (winsorize(sample_test1$`SSC.H/SSC.A`,p=c(0.005, 0.995))))
  dW1 <- density(sample_test1$`SSC.H/SSC.A`)
  plot(dW1,xlim=c(1,3.5),ylim=c(0,2.5),cex=1,
       main="Before & after doublet gating", xlab="SSC.H(size)/SSC.A(signal)",ylab="Density of events",
       lty=2,col="blue")
  # after
  sample_test2 <- subset(sample_test,sample_test$SSC.H>=0 & sample_test$SSC.A>=0)
  sample_test2$`SSC.H/SSC.A` <- sample_test2$SSC.H/sample_test2$SSC.A
  sample_test2 <- subset(sample_test2,sample_test2$`SSC.H/SSC.A` == (winsorize(sample_test2$`SSC.H/SSC.A`,p=c(0.005, 0.995))))
  dW <- density(sample_test2$`SSC.H/SSC.A`)
  lines(dW,xlim=c(1,3.5),ylim=c(0,2.5),cex=0, col="black")
  legend(3, 2.5, legend=c("before", "after"),
         col=c("blue", "black"), lty=2:1, cex=0.5,
         box.lty=0)
# doublet removal color plots before & after
  # before
  x <- densCols(sample_test0$SSC.A,sample_test0$SSC.H) 
  sample_test0$dens <- col2rgb(x)[1,] + 1L
  cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(100) #(50) #low number, hot colors #high number, cold colors
  sample_test0$col <- cols[sample_test0$dens]
  plot(sample_test0$SSC.H~sample_test0$SSC.A, data=sample_test0[order(sample_test0$dens),],
       pch=20, col=col, cex=0.5,
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),
       main="Before doublet gating",
       xlab="SSC-A (event signal)",
       ylab="SSC-H (event size)")
  # after
  x <- densCols(sample_test$SSC.A,sample_test$SSC.H) 
  sample_test$dens <- col2rgb(x)[1,] + 1L
  cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(100) #(50) #low number, hot colors #high number, cold colors
  sample_test$col <- cols[sample_test$dens]
  plot(sample_test$SSC.H~sample_test$SSC.A, data=sample_test[order(sample_test$dens),],
       pch=20, col=col, cex=0.5,
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),
       main="After doublet gating",
       xlab="SSC-A (event signal)",
       ylab="SSC-H (event size)")
  mtext(paste("Filtered sample: ", samplename0, sep=""), side = 3, outer = TRUE, cex = 1)
  
  }

graphics.off()
   
## OUTPUT TABLES POST-GATING
head(YOURDATA.CLEAN.IND)
mypath="D:/xxx/"

for (i in 1:xxx) {
  samplename = names(YOURDATA.CLEAN.IND)[i]
  sample_test <- get(samplename,YOURDATA.CLEAN.IND)
  
  write.csv(sample_test,file=paste(mypath,samplename,"_postgating",".csv",sep=""),sep=" ",col.names = TRUE,row.names = FALSE)
}

