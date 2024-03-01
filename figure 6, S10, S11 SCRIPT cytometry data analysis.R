### CODOVIREVOL - shbleSTOP project ###
### Flow cytometry analyses: FLUO ANALYSIS
### Last update by PPB June2023

###########
## SETUP ##
###########
rm(list=ls())
library('gridExtra')
library('ggplot2')
library(Rmisc)
library(RColorBrewer)
library(zoo)
library(viridis)
library(colormap)
library(plyr)
library(Hmisc)
library(lsr)
library(multcompView)
library(reshape)
library(reshape2)
library(stringr)
library(scales)
library(ggridges)
library(forcats)
library(multcomp)
#library(spatialEco)
#######################################
#####DEFINE HERE YOUR COLOR CODE######
#######################################

#for the 20 constructs
#color for the mock + empty + 13shXs
col=c("black","#666666","#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02","#bb4732","#80c861","#008080","#960018","#ffe4c4","#3f2a32","#6082b6")
#color for the mock + empty + 13 construct with replacment by mutSS
col=c("black","#666666","#1B9E77", "#D95F02", "#7570B3", "#b22d71", "#66A61E", "#ae881b","#a14433","#80c861","#008080","#712834","#ffe4c4","#3f2a32","#626e7f")
#color for the mock + empty + 18 construct 
col=c("black","#666666","#1B9E77","#D95F02","#7570B3","#E7298A","#b22d71", "#66A61E", "#E6AB02", "#ae881b","#bb4732","#a14433","#80c861","#008080","#960018","#712834","#ffe4c4","#3f2a32","#6082b6","#626e7f")



######################
## CREATE BIG TABLE ##
######################

## Define input data

#setwd("/media/marion/Commun/Documents/D2paper/Synthese_D2/R_codes_forGithub/Cytometry_Analyses/input/tables_post_gating") # The working directory should contain all the filtered CSV files (see GATING SCRIPT)
temp = list.files(pattern="*.csv")
YOURDATA.filteredinput = lapply(temp, read.csv)
names(YOURDATA.filteredinput) <-(temp)
names(YOURDATA.filteredinput)
length(YOURDATA.filteredinput) 
head(YOURDATA.filteredinput) 


#changer column names
canard<-function(i){
  colnames(i)<-c("FSC.H", "SSC.H", "GFP.H", "FSC.A", "SSC.A", "GFP.A", "Width", "Time")
  return(i)
}
YOURDATA<-lapply(YOURDATA,canard)

## Create replicates and construct columns

tag <- function(S, filename){
  S$replicate<- word(filename, 1, sep = fixed("_")) # get the first word 
  S$construct<- word(filename, 2, sep = fixed("_")) # get the second word
  S$version<- word(filename, 3, sep = fixed("_"))# get the third world 
  S$splicing<- word(filename, 4, sep = fixed("_"))
  S$mutSS<- word(filename, 5, sep = fixed("_"))
  return(S)
}

YOURDATA.tagged <- Map(tag, YOURDATA.filteredinput, names(YOURDATA.filteredinput)) 
names(YOURDATA.tagged)
length(YOURDATA.tagged) 

## Filter on number of events: discard samples with less than 30,000 events

Count.events <- function(c){
  nrow(c)}
YOURDATA.events <- lapply(YOURDATA.tagged, Count.events)
YOURDATA.lownumber <- subset(YOURDATA.events, as.numeric(YOURDATA.events)<30000)
YOURDATA.lownumber 

YOURDATA.1 <- subset(YOURDATA.tagged,names(YOURDATA.tagged)!="r1_sh2_stop_N_x_.csv"
                     &names(YOURDATA.tagged)!="r1_sh3_stop_N_x_.csv"&names(YOURDATA.tagged)!="r1_sh4_stop_Y_x_.csv"
                     &names(YOURDATA.tagged)!="r1_sh8_stop_N_x_.csv"&names(YOURDATA.tagged)!="r1_sh9_stop_N_x_.csv"
                     &names(YOURDATA.tagged)!="r7_sh4_stop_Y_x_.csv"&names(YOURDATA.tagged)!="r7_sh6_stop_Y_x_.csv"
                     &names(YOURDATA.tagged)!="r8_sh12_stop_N_x_.csv") 

length(YOURDATA.1) 

###As i Didn't removed the files with less than 30k events, I need to remane YOURDATA.tagged as YOURDATA1

YOURDATA.1 <- YOURDATA.tagged # use this if you want to bypass the subsetting of files w/ less than 30k events

## Keep 30,000 events for each remaining sample

No.neg <-  function(p){
  subset(p,p$GFP.A>=0)
}
YOURDATA.2 <- lapply(YOURDATA.1, No.neg)
Subset.events <- function(s){
  s[sample(nrow(s),30000,replace=T), ]}  #here I used both version w/ or w/o replacment accordingly to w/ or w/o the 8 files with less than 30k events
YOURDATA.2 <- lapply(YOURDATA.2, Subset.events)
length(YOURDATA.2)
YOURDATA.events <- lapply(YOURDATA.2, Count.events)
YOURDATA.events  # should be 30000 for all samples

###after running the first data manipulation and graph, only 138 .csv from the 160 original remain
# Keep only interesting columns
IncrementalTable <- c()
for (i in 1:138) {
  samplename = names(YOURDATA.2)[i]
  sample_test <- get(samplename,YOURDATA.2)
  t <- rbind(sample_test[,c("replicate","construct","version", "splicing", "mutSS", "GFP.A")])
  IncrementalTable <- rbind(IncrementalTable,t)
}
colnames(IncrementalTable) <- c("replicate","construct","version", "splicing", "mutSS", "GFP.A")
nrow(IncrementalTable) 
head(IncrementalTable)

# Print big table
mypath="C:/Users/philippe paget/Desktop/shbleSTOP_R/"
write.table(IncrementalTable,file=paste(mypath,"shbleSTOP-FL_allreplicate_30k_w_replacment_138files",".csv",sep=""),col.names = TRUE,row.names = FALSE)

## LOAD BIG TABLE - BEGIN FROM HERE WHEN ANALYSING DATA (Input bigtable provided as SUP DATA) ##

#w/replacment 138 files
bigtable <- read.table("D:/SHBLE STOP - nonSTOP - FL/shbleSTOP_R/shbleSTOP-FL_allreplicate_30k_w_replacment_138files_w_log10.csv" ,header=TRUE) 
bigtable <- shbleSTOP.FL_allreplicate_30k_w_replacment_138files_w_log10
#import verification
head(bigtable)
nrow(bigtable)
bigtable$log10_GFP.A <- log10(as.numeric(paste((bigtable$"GFP.A"))))
head(bigtable)

#I don't know why but there is infinite values in the df, though I removed neg value before doing the log transform
bigtable<- subset(bigtable, bigtable$log10_GFP.A>=0)
nrow(bigtable)

#first graph of cell pop distribution

bigtable_mock <- subset(bigtable,bigtable$construct=="mock")
p0 <- ggplot(bigtable_mock) +
  geom_density(aes(x=log10_GFP.A, color=replicate)) + 
  geom_density(aes(x=log10_GFP.A),linetype = "dashed",linewidth=0.5, color="gray", fill="black",alpha=.2) +
  ggtitle("Fluorescence signal per replicate - mock") +
  xlab("log10(GFP-A)") +
  ylab("Event density") +
  theme(plot.title= element_text(face="bold",size=16), axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(face="bold",size=14),axis.title.y = element_text(face="bold",size=14))
#geom_vline(xintercept = log10(pos_value), linetype="dashed", size=0.5)
p0

bigtable_empty <- subset(bigtable,bigtable$construct=="empty")
p01 <- ggplot(bigtable_empty) +
  geom_density(aes(x=log10_GFP.A, color=replicate)) + 
  geom_density(aes(x=log10_GFP.A),linetype = "dashed",linewidth=0.5, color="gray", fill="black",alpha=.2) +
  ggtitle("Fluorescence signal per replicate - empty") +
  xlab("log10(GFP-A)") +
  ylab("Event density") +
  theme(plot.title= element_text(face="bold",size=16), axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(face="bold",size=14),axis.title.y = element_text(face="bold",size=14))
#geom_vline(xintercept = log10(pos_value), linetype="dashed", size=0.5)
p01

#to plot every replicate from a given construct
ggplot(subset(bigtable,bigtable$construct=="sh1s")) +
  geom_density(aes(x=log10_GFP.A, color=replicate)) + 
  geom_density(aes(x=log10_GFP.A),linetype = "dashed",linewidth=0.5, color="gray", fill="black",alpha=.2) +
  ggtitle("Fluorescence signal per replicate - sh1s") +
  xlab("log10(GFP-A)") +
  ylab("Event density") +
  theme(plot.title= element_text(face="bold",size=16), axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(face="bold",size=14),axis.title.y = element_text(face="bold",size=14))
  #geom_vline(xintercept = log10(pos_value), linetype="dashed", size=0.5)


########to plot every construct from a given replicate

## Reordonner le tableau pour mettre les conditions dans l'ordre voulu
bigtable$construct <- factor(bigtable$construct, levels = c("mock","empty","sh1s",
                                                            "sh2s","sh3s","sh4s","sh4s-mutSS",
                                                            "sh5s","sh6s","sh6s-mutSS","sh7s","sh7s-mutSS","sh8s",
                                                            "sh9s","sh10s","sh10s-mutSS", "sh11s", "sh12s", "sh13s", "sh13s-mutSS"))


#to plot each replicate
ggplot(subset(bigtable,bigtable$replicate=="r1")) +
  geom_density(aes(x=log10_GFP.A, color=construct)) + 
  ggtitle("Fluorescence signal per sample - r1") +
  xlab("log10(GFP-A)") +
  ylab("Event density") +
  theme(plot.title= element_text(face="bold",size=16), axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(face="bold",size=14),axis.title.y = element_text(face="bold",size=14))

#to plot the stop version against the mutSS version
ggplot(subset(bigtable,bigtable$splicing=="Y")) +
  geom_density(aes(x=log10_GFP.A, color=mutSS)) + 
  ggtitle("Fluorescence signal stop vs. mutSS") +
  xlab("log10(GFP-A)") +
  ylab("Event density") +
  theme(plot.title= element_text(face="bold",size=16), axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(face="bold",size=14),axis.title.y = element_text(face="bold",size=14))


#to plot a given construct with or without mutSS
ggplot(subset(bigtable,bigtable$construct=="sh4s")) +
  geom_density(aes(x=log10_GFP.A, color=mutSS)) + 
  ggtitle("Fluorescence signal of sh4s stop vs. mutSS") +
  geom_vline(xintercept=c(3,4.261335), linetype="dashed")+
  xlab("log10(GFP-A)") +
  ylab("Event density") +
  theme(plot.title= element_text(face="bold",size=16), axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(face="bold",size=14),axis.title.y = element_text(face="bold",size=14))

#to plot each spliced construct with or w/o mutSS
bigtable_splice <- subset(bigtable,bigtable$splicing=="Y")
ggplot(bigtable_splice) +
  geom_density(aes(x=log10_GFP.A, color=mutSS)) + 
  ggtitle("Fluorescence signal of shble stop vs. mutSS") +
  geom_vline(xintercept=c(4.261335), linetype="dashed")+
  facet_wrap(.~construct)+
  xlab("log10(GFP-A)") +
  ylab("Event density") +
  theme(plot.title= element_text(face="bold",size=16), axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(face="bold",size=14),axis.title.y = element_text(face="bold",size=14))

#######################################################################################
######################my attempts on HPV data to generate graphs to compare with marion

#plot the average of all replicate for a given construct
ggplot(bigtable[bigtable$construct == "mock",], aes(x= log10_GFP.A, y = construct, fill= construct)) +
  geom_density_ridges() +
  theme_ridges()

#one plot for each construct with every replicate on y axis
ggplot(bigtable, aes(x= log10_GFP.A, y=replicate, color=replicate, fill=replicate)) +
  geom_density_ridges(color="black", alpha = 0.4) +
  geom_vline(xintercept=c(3,4.261335), linetype="dashed") + # change modulo threshold .99 mock
  facet_wrap(.~construct) +
  labs(   x = "log10(GFP-A)", y = "Event density") +
  theme_ridges(center_axis_labels = T)

#one plot for mutSS vs plot for spliced construct
ggplot(bigtable, aes(x= log10_GFP.A, y = replicate, color=replicate, fill=replicate)) +
  geom_density_ridges(color="black", alpha = 0.5) +
  geom_vline(xintercept=c(3,4.261335), linetype="dashed") +
  facet_wrap(.~mutSS) +
  labs(   x = "log10(GFP-A)", y = "Event density") +
  theme_ridges(center_axis_labels = T)

#one plot for sh4 mutss and another for spliced
bigtable_sh4s <- subset(bigtable,bigtable$construct=="sh4s")
bigtable_sh4s_mutSS <- subset(bigtable,bigtable$construct=="sh4s-mutSS")
bigtable_sh4 <- rbind(bigtable_sh4s, bigtable_sh4s_mutSS)


ggplot(bigtable_sh4, aes(x= log10_GFP.A, y = construct)) +
  geom_density_ridges(color="black", alpha = 0.5) +
  geom_vline(xintercept=c(3,4.261335), linetype="dashed") +
  labs(   x = "log10(GFP-A)", y = "Event density") +
  theme_ridges(center_axis_labels = T)

#color by quartile
ggplot(bigtable_sh4, aes(x = log10_GFP.A, y = construct, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles =  c(0.5, 0.8, 0.9), quantile_lines=T
  ) +
  scale_fill_manual(
    name = "Quantiles", values = c("grey", "grey", "lightgreen", "green"),
    labels = c("(0, 0.5]", "(0.5, 0.8)", "(0.8,0.9]","(0.9,1]")
  )+
  geom_vline(xintercept=c(3,4.261335), linetype="dashed") +
  labs(   x = "log10(GFP-A)", y = "Event density") +
  theme_ridges(center_axis_labels = T)

#for every construct doing splicing
bigtable_splicing <- subset(bigtable,bigtable$splicing=="Y")
bigtable_splicing$construct <- factor(bigtable_splicing$construct, levels = c("sh4s","sh4s-mutSS",
                                                          "sh6s","sh6s-mutSS",
                                                          "sh7s","sh7s-mutSS",
                                                            "sh10s","sh10s-mutSS",
                                                            "sh13s", "sh13s-mutSS"))

ggplot(bigtable_splicing, aes(x = log10_GFP.A, y = construct, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles =  c(0.5, 0.8, 0.9), quantile_lines=T
  ) +
  scale_fill_manual(
    name = "Quantiles", values = c("grey", "grey", "lightgreen", "green"),
    labels = c("(0, 0.5]", "(0.5, 0.8)", "(0.8,0.9]","(0.9,1]")
  )+
  geom_vline(xintercept=c(3,4.261335), linetype="dashed") +
  labs(   x = "log10(GFP-A)", y = "Event density") +
  theme_ridges(center_axis_labels = T)

###exemple de grepl if needed
#ggplot(bigtable[grepl("16", bigtable$condition),], aes(x= log10_GFP.A, y = replicate, color=replicate, fill=replicate)) +
 # geom_density_ridges(color="black",alpha = 0.5) +
 # geom_vline(xintercept=c(3,4.575084), linetype="dashed") +
 #facet_wrap(.~construct) +
 # labs(   x = "log10(GFP-A)", y = "cell count") +
  #theme_ridges(center_axis_labels = T)



#############################
# calculation of the value of "positivity" - based on all events
bigtable_mock <-subset(bigtable,bigtable$construct=="mock")
density(bigtable_mock$GFP.A)
pos_value <- quantile(bigtable_mock$GFP.A, probs = seq(0.99,1))
pos_value #for r1,r2,r3,r5,r6 and r7 mock = 18253.04
pos_value_log10<-quantile(bigtable_mock$log10_GFP.A, probs = seq(0.99,1))
pos_value_log10 #4.261335 --> included in prior graphs 


##########################
#at this step wipe everything and load the file again
bigtable <- read.table("C:/Users/philippe paget/Desktop/shbleSTOP_R/shbleSTOP-FL_allreplicate_30k_w_replacment_138files_w_log10.csv" ,header=TRUE) 
#added the log transform which was not saved
mypath="C:/Users/philippe paget/Desktop/shbleSTOP_R/"
write.table(bigtable,file=paste(mypath,"shbleSTOP-FL_allreplicate_30k_w_replacment_138files_w_log10",".csv",sep=""),col.names = TRUE,row.names = FALSE)


#I kept marion's nomenclature to be sure (for her noout was her original bigtable minus to bad replicates but I already discard what I don't want)
bigtable_noout <-bigtable

# Create tables per construct
bigtable_noout_sh1s <- subset(bigtable_noout,bigtable_noout$construct=="sh1s")
bigtable_noout_sh2s <- subset(bigtable_noout,bigtable_noout$construct=="sh2s")
bigtable_noout_sh3s <- subset(bigtable_noout,bigtable_noout$construct=="sh3s")
bigtable_noout_sh4s <- subset(bigtable_noout,bigtable_noout$construct=="sh4s")
bigtable_noout_sh4s_mutSS <- subset(bigtable_noout,bigtable_noout$construct=="sh4s-mutSS")
bigtable_noout_sh5s <- subset(bigtable_noout,bigtable_noout$construct=="sh5s")
bigtable_noout_sh6s <- subset(bigtable_noout,bigtable_noout$construct=="sh6s")
bigtable_noout_sh6s_mutSS <- subset(bigtable_noout,bigtable_noout$construct=="sh6s-mutSS")
bigtable_noout_sh7s <- subset(bigtable_noout,bigtable_noout$construct=="sh7s")
bigtable_noout_sh7s_mutSS <- subset(bigtable_noout,bigtable_noout$construct=="sh7s-mutSS")
bigtable_noout_sh8s <- subset(bigtable_noout,bigtable_noout$construct=="sh8s")
bigtable_noout_sh9s <- subset(bigtable_noout,bigtable_noout$construct=="sh9s")
bigtable_noout_sh10s <- subset(bigtable_noout,bigtable_noout$construct=="sh10s")
bigtable_noout_sh10s_mutSS <- subset(bigtable_noout,bigtable_noout$construct=="sh10s-mutSS")
bigtable_noout_sh11s <- subset(bigtable_noout,bigtable_noout$construct=="sh11s")
bigtable_noout_sh12s <- subset(bigtable_noout,bigtable_noout$construct=="sh12s")
bigtable_noout_sh13s <- subset(bigtable_noout,bigtable_noout$construct=="sh13s")
bigtable_noout_sh13s_mutSS <- subset(bigtable_noout,bigtable_noout$construct=="sh13s-mutSS")
bigtable_noout_mock <- subset(bigtable_noout,bigtable_noout$construct=="mock")
bigtable_noout_empty <- subset(bigtable_noout,bigtable_noout$construct=="empty")

bigtable_noout$construct <- factor(bigtable_noout$construct, levels = c("mock","empty","sh1s",
                                                            "sh2s","sh3s","sh4s","sh4s-mutSS",
                                                            "sh5s","sh6s","sh6s-mutSS","sh7s","sh7s-mutSS","sh8s",
                                                            "sh9s","sh10s","sh10s-mutSS", "sh11s", "sh12s", "sh13s", "sh13s-mutSS"))


                     #########################################
                     ## DISTRIBUTION OF FLUORESCENCE SIGNAL ##
                     #########################################

p_all_noout <- ggplot(bigtable_noout, aes(x=log10_GFP.A, color=construct, fill=construct)) +
  geom_density(alpha=.2) +
  xlab("Fluorescence level [log10(FITC-A)]") +
  ylab("Event density") +
  scale_fill_manual(values=col,name = "sh gene versions", labels = c("Mock","Superempty","Empty","1: Common", "2: Common - GC rich", "3: Common - AT rich","4: Rare", "5: Rare - GC rich", "6: Rare - AT rich"))+
  scale_color_manual(values=col,name = "sh gene versions", labels = c("Mock","Superempty","Empty","1: Common", "2: Common - GC rich", "3: Common - AT rich","4: Rare", "5: Rare - GC rich", "6: Rare - AT rich"))+
  theme_bw()+
  theme(legend.position=c(0.85,0.6), axis.text.x = element_text(size=14),axis.text.y = element_text(size=14),axis.title.x = element_text(face="bold",size=16),axis.title.y = element_text(face="bold",size=16))+
  geom_vline(xintercept = log10(pos_value), linetype="dashed", linewidth=0.5)
p_all_noout 

## all events - per line
#library(forcats)
bigtable_noout$construct=fct_relevel(bigtable_noout$construct,c("mock","empty","sh6","sh5","sh4","sh3","sh2","sh1"))

library(ggridges)
head(bigtable_noout)
ggplot(bigtable_noout, aes(x = log10_GFP.A, y = construct, fill = construct,height = ..density..)) +
  geom_density_ridges(na.rm = FALSE,inherit.aes = TRUE,panel_scaling = FALSE,mapping = NULL,data = NULL,stat = "density")+ 
  theme_ridges() + 
  xlab("Fluorescence level [log10(FITC-A)]")+
  scale_fill_manual(values=col_ridge,name = "sh gene versions", labels = c("Mock","Superempty","Empty","1: Common", "2: Common - GC rich", "3: Common - AT rich","4: Rare", "5: Rare - GC rich", "6: Rare - AT rich"))+
  scale_color_manual(values=col_ridge,name = "sh gene versions", labels = c("Mock","Superempty","Empty","1: Common", "2: Common - GC rich", "3: Common - AT rich","4: Rare", "5: Rare - GC rich", "6: Rare - AT rich"))+
  theme(legend.position = "none")+
  geom_vline(xintercept = log10(pos_value), linetype="dashed", size=0.5)

## per construct
p_all<- ggplot(bigtable_noout) +
  geom_density(aes(x=log10_GFP.A, color=replicate)) + 
  geom_density(aes(x=log10_GFP.A),linetype = "dashed",linewidth=0.5, fill="black",alpha=.2) +
  xlab("Fluorescence level [log10(FITC-A)]") +
  ylab("Event density") +
  #theme(legend.position = "none")  +
  facet_wrap(c(~construct), ncol=5, nrow=4)
  p_all

## statistical tests
sh1s=sample((subset(bigtable_noout,bigtable_noout$construct=="sh1s"&bigtable_noout$log10_GFP.A!="-Inf"))$log10_GFP.A,100000)
#I don't know why bhy even after removing neg values two time they still pop in my sampling
#to conserve the same number of event in each sampling, I choose to remove AGAIN the neg values from sh6s-mutSS but prior to sampling
#that is why she was using the subset(df, df$col==!"construct")
sh1s <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh1s"))$log10_GFP.A,100000)
sh2s <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh2s"))$log10_GFP.A,100000)
sh3s <-sample((subset(bigtable_noout,bigtable_noout$construct=="sh3s"))$log10_GFP.A,100000)
sh4s <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh4s"))$log10_GFP.A,100000)
sh4s_mutSS <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh4s-mutSS"))$log10_GFP.A,100000)
sh5s <-sample((subset(bigtable_noout,bigtable_noout$construct=="sh5s"))$log10_GFP.A,100000)
sh6s <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh6s"))$log10_GFP.A,100000)
###############################
sh6s_mutSS <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh6s-mutSS"))$log10_GFP.A,100000)
#this one poses a problem so by giving one -Inf value to remove before sampling
nrow(bigtable_noout_sh6s_mutSS) #120k (30k * 4 replicates)
bigtable_noneg_sh6s_mutSS <- subset(bigtable_noout_sh6s_mutSS,bigtable_noout_sh6s_mutSS$log10_GFP.A>=0)
nrow(bigtable_noneg_sh6s_mutSS) #119999 
#then sampling
sh6s_mutSS <- sample((subset(bigtable_noneg_sh6s_mutSS,bigtable_noneg_sh6s_mutSS$construct=="sh6s-mutSS"))$log10_GFP.A,100000)
###############################
sh7s <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh7s"))$log10_GFP.A,100000)
sh7s_mutSS <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh7s-mutSS"))$log10_GFP.A,100000)
sh8s <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh8s"))$log10_GFP.A,100000)
sh9s <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh9s"))$log10_GFP.A,100000)
sh10s <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh10s"))$log10_GFP.A,100000)
sh10s_mutSS <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh10s-mutSS"))$log10_GFP.A,100000)
sh11s <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh11s"))$log10_GFP.A,100000)
sh12s <- sample((subset(bigtable_noout,bigtable_noout$construct=="sh12s"))$log10_GFP.A,100000)
sh13s <-sample((subset(bigtable_noout,bigtable_noout$construct=="sh13s"))$log10_GFP.A,100000)
sh13s_mutSS <-sample((subset(bigtable_noout,bigtable_noout$construct=="sh13s-mutSS"))$log10_GFP.A,100000)
mock <- sample((subset(bigtable_noout,bigtable_noout$construct=="mock"))$log10_GFP.A,100000)
empty <- sample((subset(bigtable_noout,bigtable_noout$construct=="empty"))$log10_GFP.A,100000)

library(kSamples)
#distribution normality test (N must be inf to 5000)
shapiro.test(sample((bigtable_noout_mock)$log10_GFP.A, 3000)) #not normal
shapiro.test(sample((bigtable_noout_empty)$log10_GFP.A, 3000)) #not normal

# non-parametric Kruskal & Wallis
qn.test(sh1s,sh2s,sh3s,sh4s,sh5s,sh6s,sh7s,sh8s,sh9s,sh10s,sh11s,sh12s,sh13s, data = NULL, test = c("KW"),
          +       method = c("asymptotic"), dist = FALSE, Nsim = 10000)

qn.test(sh1s,sh2s,sh3s,sh4s_mutSS,sh5s,sh6s_mutSS,sh7s_mutSS,sh8s,sh9s,sh10s_mutSS,sh11s,sh12s,sh13s_mutSS, data = NULL, test = c("KW"),
          +       method = c("asymptotic"), dist = FALSE, Nsim = 10000)

#mock vs. construct
ad.test(sh1s,mock, data = NULL, method = c("asymptotic"), dist = FALSE, Nsim = 10000)
ad.test(sh2s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh3s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh4s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh4s_mutSS, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh5s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh6s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh6s_mutSS, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh7s,mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh7s_mutSS, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh8s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh9s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh10s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh10s_mutSS, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh11s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh12s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh13s, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh13s_mutSS, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(mock, mock, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(mock, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)

#empty vs. construct
ad.test(sh1s,empty, data = NULL, method = c("asymptotic"), dist = FALSE, Nsim = 10000)
ad.test(sh2s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh3s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh4s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh4s_mutSS, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh5s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh6s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh6s_mutSS, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh7s,empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh7s_mutSS, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh8s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh9s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh10s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh10s_mutSS, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh11s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh12s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh13s, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh13s_mutSS, empty, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)

#spliced vs. mutSS
ad.test(sh4s, sh4s_mutSS, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh6s, sh6s_mutSS, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh7s, sh7s_mutSS, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh10s, sh10s_mutSS, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)
ad.test(sh13s, sh13s_mutSS, data = NULL, method = c("asymptotic"),dist = FALSE, Nsim = 10000)


                          ###################  
                          ## GAUSSIAN FITS ##  
                          ###################

bigtable <- read.table("C:/Users/philippe paget/Desktop/shbleSTOP_R/shbleSTOP-FL_allreplicate_30k_w_replacment_138files_w_log10.csv" ,header=TRUE) 
bigtable_noout <- bigtable

## LOAD BIG TABLE - BEGIN FROM HERE WHEN ANALYSING DATA (Input bigtable provided as SUP DATA) ##

#w/replacment 138 files
bigtable <- read.table("D:/SHBLE STOP - nonSTOP - FL/shbleSTOP_R/shbleSTOP-FL_allreplicate_30k_w_replacment_138files_w_log10.csv" ,header=TRUE) 
bigtable <- shbleSTOP.FL_allreplicate_30k_w_replacment_138files_w_log10
#import verification
head(bigtable)
nrow(bigtable)
bigtable$log10_GFP.A <- log10(as.numeric(paste((bigtable$"GFP.A"))))
head(bigtable)

#I don't know why but there is infinite values in the df, though I removed neg value before doing the log transform
bigtable<- subset(bigtable, bigtable$log10_GFP.A>=0)
nrow(bigtable)

bigtable_noout_empty <- subset(bigtable, bigtable$construct=="empty")
bigtable_noout_sh4s <- subset(bigtable, bigtable$construct=="sh4s")
## MCLUST graphs
set.seed(1) 
sample_empty<-sample(bigtable_noout_empty$log10_GFP.A,30000)
sample_sh4s<-sample(bigtable_noout_sh4s$log10_GFP.A,30000)
sample_mock<-sample(bigtable_noout_mock$log10_GFP.A,30000)

library(mclust)
set.seed(1)
test <-densityMclust(sample_sh4s,modelNames = "V") # V = variable/unqual variance (one-dimensional)
plot(test, what="BIC")

plotDensityMclust1(test,data=sample_sh4s,breaks=100)
test$parameters
test$parameters$pro
test$parameters$mean
summary(test, parameters=TRUE)
densityMclust.diagnostic(test, type = c("cdf", "qq"), col=c("blue","red"),lty=c(3,1) , legend=T) ###supdata


## MCLUST tables HERE YOU DEFINE THE SAMPLES YOU WANT TO INCLUDE IN THE TEST
head(bigtable_noout)
nrow(bigtable_noout)
table_formclust <- subset(bigtable_noout,bigtable_noout$log10_GFP.A!=-Inf & bigtable_noout$construct!="mock")
table_formclust <- subset(bigtable_noout,bigtable_noout$log10_GFP.A!=-Inf & bigtable_noout$construct!="mock" & bigtable_noout$construct!="empty")
table_formclust <- subset(bigtable_noout,bigtable_noout$log10_GFP.A!=-Inf & bigtable_noout$construct!="mock" & bigtable_noout$construct!="empty"
                         & bigtable_noout$construct!="sh4s_mutSS" & bigtable_noout$construct!="sh6s_mutSS" & bigtable_noout$construct!="sh7s_mutSS"
                         & bigtable_noout$construct!="sh10s_mutSS" & bigtable_noout$construct!="sh13s_mutSS")
table_formclust <- subset(bigtable_noout,bigtable_noout$log10_GFP.A!=-Inf & bigtable_noout$construct!="mock" & bigtable_noout$construct!="empty"
                          & bigtable_noout$construct!="sh4s" & bigtable_noout$construct!="sh6s" & bigtable_noout$construct!="sh7s"
                          & bigtable_noout$construct!="sh10s" & bigtable_noout$construct!="sh13s")
table_formclust <- subset(bigtable_noout,bigtable_noout$log10_GFP.A!=-Inf)

nrow(table_formclust)

mclust_table=c()
for(x in unique(bigtable$construct)){   #table_formclust
  y = subset(bigtable,bigtable$construct==x)
 for (i in unique(y$replicate)){
    j = subset(y,y$replicate==i)
    r = as.numeric(paste(j$log10_GFP.A))
    set.seed(1)
    k=densityMclust(r,G=2,na.rm = TRUE,modelNames = "V")
    mean1 = k$parameters$mean[1]
    mean2 = k$parameters$mean[2]
    var1 = k$parameters$variance$sigmasq[1]
    var2 = k$parameters$variance$sigmasq[2]
    tmp = data.frame(construct=x, replicate = i,mean1=mean1, mean2=mean2, var1=var1, var2=var2)
    mclust_table = rbind(mclust_table, tmp)
  }
}
head(mclust_table)
nrow(mclust_table)
mypath= "E:/article SHBLEstop 240124 1814/"   #"C:/Users/philippe paget/Desktop/shbleSTOP_R/"
write.table(mclust_table,file=paste(mypath,"mclust_mean_var",".csv",sep=""),col.names = TRUE,row.names = FALSE)

# median of the means per construct
mean_table=c()
for(x in unique(mclust_table$construct)){
  y = subset(mclust_table,mclust_table$construct==x)
  mean1_mean = median(y$mean1)
  mean2_mean = median(y$mean2)
  tmp = data.frame(construct=x, median_mean1=mean1_mean, median_mean2=mean2_mean)
  mean_table = rbind(mean_table, tmp)
}
mean_table
write.table(mean_table,file=paste(mypath,"mclust_median_ofthemean_ofeachreplicate_20constructs",".csv",sep=""),col.names = TRUE,row.names = FALSE)

#the result matrix generated by the next loop don't accept "-" in the row names so:
mclust_table$construct <- gsub('-', '', mclust_table$construct)


# mean 1
statTest_mean1 = c()
for(i in unique(mclust_table$construct)){
  p = pairwise.wilcox.test(as.numeric(paste(mclust_table$mean1)), mclust_table$construct, p.adjust.method = "none") # test is here
  l = melt(p$p.value)[,3]
  print(l) # print the matrix of p-values
  names(l) = paste(melt(p$p.value)[,2], melt(p$p.value)[,1], sep="-")
  pvalues_mean1 <- data.frame(l) #added this to retrieve the p-values
  l = l[which(is.na(l)==FALSE)]
  t = multcompLetters(l,  compare='<', threshold = 0.05, Letters = c(letters, LETTERS, "."))
  tmp = data.frame(
    construct = names(t$Letters),
    letter = as.character(t$Letters),
    c = rep(i, length(t$Letters))
  )
  statTest_mean1 = rbind(statTest_mean1, tmp)
}
statTest_mean1$construct <- factor(statTest_mean1$construct, c("empty","sh1s",
                                                               "sh2s","sh3s","sh4s","sh4smutSS",
                                                               "sh5s","sh6s","sh6smutSS","sh7s","sh7smutSS","sh8s",
                                                               "sh9s","sh10s","sh10smutSS", "sh11s", "sh12s", "sh13s", "sh13smutSS"))
write.table(statTest_mean1,file=paste(mypath,"pairwiseWilcoxtest_nocorec_mean1_20constructs",".csv",sep=""),col.names = TRUE,row.names = FALSE)
write.table(pvalues_mean1,file=paste(mypath,"pairwiseWilcoxtest_nocorec_mean1_20constructs_pvalues",".csv",sep=""),col.names = TRUE,row.names = TRUE)


              ######################################################################
              ###                   double box-plot Fig 6                          #
              ######################################################################

#load your library of choice
library(gridExtra)
library(ggplot2)
library(ggExtra)
library(Rmisc)
library(RColorBrewer)
library(zoo)
library(viridis)
library(colormap)
library(plyr)
library(Hmisc)
library(lsr)
library(multcompView)
library(reshape)
library(reshape2)
library(stringr)
library(scales)
library(ggridges)
library(forcats)
library(multcomp)
library(janitor)
library(hrbrthemes)
library(ggpubr)
#####################some preanalyses basic stuff

#import your df with push-buton system (avoid to ruin your header/ special characters during importation)

#DONT FORGET TO TRANSFORM YOUR DATA COL INTO NUM OR INTEGER !!!!!!
big_recap <- mclust_mean_var_PCAaxis
big_recap[,9:18] <- lapply(big_recap[,9:18], as.numeric)

#typically, we don't include mock and empty conditions in our graph
big_recap <- subset(big_recap, big_recap$construct!="mock")
big_recap <- subset(big_recap, big_recap$construct!="empty")


###
# e) choose the label and the order you want your conditions to appear on your graph, label and order (factor) have to match if you want correct graph plotting !, choose a color palette for your graphs
#13shble 
big_recap <- subset(big_recap, big_recap$cub=="yes")
big_recap$construct<- factor(big_recap$construct, c("sh1","sh2","sh3","sh4","sh5","sh6","sh7","sh8","sh9","sh10","sh11","sh12","sh13"))
label <- c("sh1","sh2","sh3","sh4","sh5","sh6","sh7","sh8","sh9","sh10","sh11","sh12","sh13")
#        1         2           3         4        5         6       7           8         9        10         11        12      13  
col=c("#D95F02","#1B9E77","#6082b6", "#ae881b","#909495","#a14433","#66A61E","#332298","#CB2030","#b22d71","#3f2a32","#d1ab60","#7570B3")

#10 splice
big_recap <- subset(big_recap, big_recap$splicing=="yes")
big_recap$construct <- factor(big_recap$construct, c("sh4sp","sh6sp","sh7sp","sh10sp","sh13sp","sh4","sh6","sh7","sh10","sh13"))
label = c("sh4sp","sh6sp","sh7sp","sh10sp","sh13sp","sh4","sh6","sh7","sh10","sh13")
col=c("#E6AB02", "#ae881b","#bb4732","#a14433", "#80c861","#66A61E","#E7298A","#b22d71", "#6082b6","#7570B3")

#####FIRST STAT TEST
statTest_X1= c()
for(i in unique(big_recap$construct)){
  p = pairwise.wilcox.test(as.numeric(paste(big_recap$mean1)), big_recap$construct,paired=F, p.adjust.method = "BH") # test is here
  l = melt(p$p.value)[,3]
  print(l)
  names(l) = paste(melt(p$p.value)[,2], melt(p$p.value)[,1], sep="-")
  l = l[which(is.na(l)==F)]
  pvalues_stat_test_X_woNA <- data.frame(l) #added this to retreive the pvalues w/o NaN
  t = multcompLetters(l,  compare='<', threshold = 0.05, Letters = c(letters, LETTERS, "."))
  tmp = data.frame(
    construct = names(t$Letters),
    letter = as.character(t$Letters),
    c = rep(i, length(t$Letters))
  )
  statTest_X1 = rbind(statTest_X1, tmp)
}

statTest_X2 = c()
for(i in unique(big_recap$construct)){
  p = pairwise.wilcox.test(as.numeric(paste(big_recap$mean2)), big_recap$construct,paired=F, p.adjust.method = "BH") # test is here
  l = melt(p$p.value)[,3]
  print(l)
  names(l) = paste(melt(p$p.value)[,2], melt(p$p.value)[,1], sep="-")
  l = l[which(is.na(l)==F)]
  pvalues_stat_test_X_woNA <- data.frame(l) #added this to retreive the pvalues w/o NaN
  t = multcompLetters(l,  compare='<', threshold = 0.05, Letters = c(letters, LETTERS, "."))
  tmp = data.frame(
    construct = names(t$Letters),
    letter = as.character(t$Letters),
    c = rep(i, length(t$Letters))
  )
  statTest_X2 = rbind(statTest_X2, tmp)
}


statTest_X1$construct <- factor(statTest_X1$construct, label) #the A to Z filtering always put 10 before 1 so you need to reorder the construct to start w/ 1,2,3....viva l'alg?rie
statTest_X1
statTest_X2$construct <- factor(statTest_X2$construct, label) #the A to Z filtering always put 10 before 1 so you need to reorder the construct to start w/ 1,2,3....viva l'alg?rie
statTest_X2

#############The FIRST graph
tiff("SP - mean of gaussian double box plot - wo stat.tif", units="in", width=8, height=6, res=300, compression = 'none')

col=c("#63e263","#1B911C")

p_X <- ggplot(big_recap, aes(x=construct, y= gaussian, fill = gaussian_num))+theme_ipsum()

p_X + geom_boxplot(width=0.8)+
  theme(legend.position="none")+
  #geom_jitter(width = 0.1, size=3)+
  ylim(min(big_recap$gaussian)-0.1,max(big_recap$gaussian)+0.1)+
  #scale_y_continuous(trans = "log10")+     #transform les data en log10 mais uniquemment pour la rpz graphique
  scale_x_discrete(labels=label)+
  scale_fill_manual(values=col,name = "shble version", labels=label)+
  scale_color_manual(values=col,name = "shble version" , labels=label)+
  theme(axis.text.x = element_text(size=16, face ="bold", color = "black", angle=45),
        axis.text.y = element_text(size=16, face ="bold", color = "black"),
        axis.title.x = element_text(face="bold",size=18, hjust=0.5),
        axis.title.y = element_text(face="bold",size=18, hjust=0.5))+
  #ggtitle("") +
  xlab("") + 
  ylab("Mean of gaussian")
geom_text(data = statTest_X1, aes(x=construct, y=max(big_recap$gaussian)+1,label=statTest_X1$letter),size=6) # in aes, add angle=90 tu turn the letters when too many

dev.off()


                #############################
                ## MEDIAN AND TOTAL SIGNAL ##
                #############################

bigtable <- read.table("D:/SHBLE STOP - nonSTOP - FL/shbleSTOP_R/shbleSTOP-FL_allreplicate_30k_w_replacment_138files_w_log10.csv" ,header=TRUE) 
bigtable_noout <- bigtable

bigtable_noout$construct <- factor(bigtable_noout$construct, c("mock","empty","sh1s",
                                   "sh2s","sh3s","sh4s","sh4s-mutSS",
                                   "sh5s","sh6s","sh6s-mutSS","sh7s","sh7s-mutSS","sh8s",
                                   "sh9s","sh10s","sh10s-mutSS", "sh11s", "sh12s", "sh13s", "sh13s-mutSS"))


## calculations
#the bigtable_noout table possess each replicate of each constrcut with 30k events
head(bigtable_noout)
stat_per_rep <- ddply(bigtable_noout, c("construct", "replicate"), summarise, sum_fluo=sum(as.numeric(paste(GFP.A)),na=TRUE), median_fluo=median(as.numeric(paste(GFP.A)),na=TRUE))
head(stat_per_rep) 

mypath="XXX"
write.table(stat_per_rep,file=paste(mypath,"shbleSTOP-FL_allreplicate_30k_w_replacment_138files_fluoSUMandMED_allreplicate",".csv",sep=""),col.names = TRUE,row.names = FALSE)


median_table=c()
for(x in unique(stat_per_rep$construct)){
  y = subset(stat_per_rep,stat_per_rep$construct==x)
  sum_fluo = median(y$sum_fluo)
  median_fluo = median(y$median_fluo)
  tmp = data.frame(construct=x, sum_fluo=sum_fluo, median_fluo=median_fluo)
  median_table = rbind(median_table, tmp)
}
median_table
write.table(median_table,file=paste(mypath,"shbleSTOP-FL_allreplicate_30k_w_replacment_138files_MEDAINof_fluoSUMandMED_fromeachreplicate",".csv",sep=""),col.names = TRUE,row.names = FALSE)


                  ###############
                  ## INTERVALS ##
                  ###############

bigtable <- read.table("D:/SHBLE STOP - nonSTOP - FL/shbleSTOP_R/shbleSTOP-FL_allreplicate_30k_w_replacment_138files_w_log10.csv" ,header=TRUE) 
bigtable_noout <- bigtable

## How many events per interval?
head(bigtable_noout)
nrow(bigtable_noout) 

#####don't generate false resultst for wilcoxtest because of the "-"
bigtable_noout$construct <- gsub('-', '', bigtable_noout$construct)

density(bigtable_noout_mock$GFP.A)
min(bigtable_noout$GFP.A)
pos_value <- quantile(bigtable_noout_mock$GFP.A, probs = seq(0.99,1))
pos_value #18253.04 
pos_value_log10<-quantile(bigtable_noout_mock$log10_GFP.A, probs = seq(0.99,1))
pos_value_log10 #4.261335

pos_table_per_replicate=c()
for(x in unique(bigtable_noout$construct)){
  y = subset(bigtable_noout,bigtable_noout$construct==x)
  for (i in unique(y$replicate)){
    j = subset(y,y$replicate==i)
    k = subset(j,j$GFP.A >= 18253.04)
    l = nrow(k)
    m = l/nrow(j)
    tmp = data.frame(construct=x, replicate = i, interval="pos", events = l,  prop = m)
    tmp
    pos_table_per_replicate = rbind(pos_table_per_replicate, tmp)
  }
}
pos_table_per_replicate

neg_table_per_replicate=c()
for(x in unique(bigtable_noout$construct)){
  y = subset(bigtable_noout,bigtable_noout$construct==x)
  for (i in unique(y$replicate)){
    j = subset(y,y$replicate==i)
    k = subset(j,j$GFP.A < 18253.04)
    l = nrow(k)
    m = l/nrow(j)
    tmp = data.frame(construct=x, replicate = i, interval="neg", events = l,  prop = m)
    tmp
    neg_table_per_replicate = rbind(neg_table_per_replicate, tmp)
  }
}
neg_table_per_replicate

i1_table_per_replicate=c()
for(x in unique(bigtable_noout$construct)){
  y = subset(bigtable_noout,bigtable_noout$construct==x)
  for (i in unique(y$replicate)){
    j = subset(y,y$replicate==i)
    k = subset(j,log10(j$GFP.A) >= 4.261335 & log10(j$GFP.A) < 5.261335)
    l = nrow(k)
    m = l/nrow(j)
    tmp = data.frame(construct=x, replicate = i, interval="i1", events = l,  prop = m)
    tmp
    i1_table_per_replicate = rbind(i1_table_per_replicate, tmp)
  }
}
head(i1_table_per_replicate)

i2_table_per_replicate=c()
for(x in unique(bigtable_noout$construct)){
  y = subset(bigtable_noout,bigtable_noout$construct==x)
  for (i in unique(y$replicate)){
    j = subset(y,y$replicate==i)
    k = subset(j,log10(j$GFP.A) >= 5.261335 & log10(j$GFP.A) < 6.261335)
    l = nrow(k)
    m = l/nrow(j)
    tmp = data.frame(construct=x, replicate = i, interval="i2", events = l,  prop = m)
    tmp
    i2_table_per_replicate = rbind(i2_table_per_replicate, tmp)
  }
}
head(i2_table_per_replicate)

i3_table_per_replicate=c()
for(x in unique(bigtable_noout$construct)){
  y = subset(bigtable_noout,bigtable_noout$construct==x)
  for (i in unique(y$replicate)){
    j = subset(y,y$replicate==i)
    k = subset(j,log10(j$GFP.A) >= 6.261335 )
    l = nrow(k)
    m = l/nrow(j)
    tmp = data.frame(construct=x, replicate = i, interval="i3", events = l,  prop = m)
    tmp
    i3_table_per_replicate = rbind(i3_table_per_replicate, tmp)
  }
}
head(i3_table_per_replicate)

#i4_table_per_replicate=c()
#for(x in unique(bigtable_noout$construct)){
#  y = subset(bigtable_noout,bigtable_noout$construct==x)
#  for (i in unique(y$replicate)){
#    j = subset(y,y$replicate==i)
#    k = subset(j,log10(j$GFP.A) >= 7)
#
#   l = nrow(k)
#    m = l/nrow(j)
#    tmp = data.frame(construct=x, replicate = i, interval="i4", events = l,  prop = m)
#    tmp
#i4_table_per_replicate = rbind(i4_table_per_replicate, tmp)
#  }
#}
#head(i4_table_per_replicate)

events_per_intervals_per_replicate=rbind(pos_table_per_replicate,neg_table_per_replicate,i1_table_per_replicate,i2_table_per_replicate,i3_table_per_replicate) #,i4_table_per_replicate)
nrow(events_per_intervals_per_replicate) #822

mypath="D:/SHBLE STOP - nonSTOP - FL/shbleSTOP_R/"
write.table(events_per_intervals_per_replicate,file=paste(mypath,"FACS_intervals_all_samples",".csv",sep=""),col.names = TRUE,row.names = FALSE)

prop_intervals_per_replicate_forgraph <- events_per_intervals_per_replicate

## giant table with all proportions
prop_intervals_per_replicate_forgraph <- read.table("D:/SHBLE STOP - nonSTOP - FL/shbleSTOP_R/FACS_intervals_all_samples.csv",header=TRUE)
head(prop_intervals_per_replicate_forgraph)

summary_table=c()
for(x in unique(prop_intervals_per_replicate_forgraph$construct)){
  y = subset(prop_intervals_per_replicate_forgraph,prop_intervals_per_replicate_forgraph$construct==x)
  for (i in unique(y$interval)){
    j = subset(y,y$interval==i)
    mean = mean(j$prop)
    sd = sd(j$prop)
    tmp = data.frame(construct=x, interval = i, sd=sd,  mean=mean)
    tmp
    summary_table = rbind(summary_table, tmp)
  }
}
summary_table


write.table(summary_table,file=paste(mypath,"FACS_intervals_per_construct_w_S.D.",".csv",sep=""),col.names = TRUE,row.names = FALSE)

################################ plot proportions - All conditions

library(forcats)
summary_table=subset(summary_table,summary_table$interval!="pos")
summary_table$interval=fct_relevel(summary_table$interval,c("neg","i1","i2","i3"))
summary_table$construct=fct_relevel(summary_table$construct,c("mock","empty","sh1s",
                                                              "sh2s","sh3s","sh4s","sh4smutSS",
                                                              "sh5s","sh6s","sh6smutSS","sh7s","sh7smutSS","sh8s",
                                                              "sh9s","sh10s","sh10smutSS", "sh11s", "sh12s", "sh13s", "sh13smutSS"))
col=c("#666666","#1B9E77","#80c861","#66A61E")


p = ggplot(data=summary_table, aes(x=interval, y=mean, fill=interval))
p + geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  geom_text(aes(label=round(mean,digit=2), y=(mean+0.4)), color="black", size=4)+ #display proportion of interval on top of hitsogram
  facet_grid(c(~construct)) +
  labs(x = "fluorescence intervals", y = "proportion of events") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(face="bold",size=14),
        strip.text.x = element_text(size = 12, colour = "black", margin = margin(.1, .1, .1, 0, "cm")),
        strip.text.y = element_text(size = 12),
        strip.background = element_rect(size=0.3))+
   scale_fill_viridis_d(option = "D")+
  scale_color_viridis_d(option = "D")
geom_text(data = statTest_neg, aes(x=construct, y=max(log10(summary_table$mean))+0.2, label=statTest_neg$letter),size=4)

# plot proportions - splice vs. mutSS
summary_table=subset(summary_table,summary_table$construct!="mock"&summary_table$construct!="empty"&summary_table$interval!="pos" &summary_table$construct!="sh1s"
                     &summary_table$construct!="sh2s"&summary_table$construct!="sh3s"&summary_table$construct!="sh5s"&summary_table$construct!="sh8s"&summary_table$construct!="sh9s"
                     &summary_table$construct!="sh11s"&summary_table$construct!="sh12s")
summary_table$interval=fct_relevel(summary_table$interval,c("neg","i1","i2","i3"))
summary_table$construct=fct_relevel(summary_table$construct,c("sh4s","sh4s-mutSS",
                                                              "sh6s","sh6s-mutSS","sh7s","sh7s-mutSS",
                                                              "sh10s","sh10s-mutSS","sh13s","sh13s-mutSS"))

p = ggplot(data=summary_table, aes(x=interval, y=mean, fill=interval))
p + geom_bar(stat="identity")+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
    geom_text(aes(label=round(mean,digit=2), y=(mean+0.4)), color="black", size=4) +
    facet_grid(c(~construct)) +
    labs(x = "fluorescence intervals", y = "proportion of events") +
    theme(axis.text=element_text(size=12), 
          axis.title=element_text(face="bold",size=14),
          strip.text.x = element_text(size = 12, colour = "black", margin = margin(.1, .1, .1, 0, "cm")),
          strip.text.y = element_text(size = 14),
          strip.background = element_rect(size=0.3))+
    scale_fill_viridis_d(option = "D")+
    scale_color_viridis_d(option = "D")

# stat per interval
neg_table_per_replicate_forstat=subset(neg_table_per_replicate,neg_table_per_replicate$interval!="pos")
i1_table_per_replicate_forstat=subset(i1_table_per_replicate,i1_table_per_replicate$interval!="pos")
i2_table_per_replicate_forstat=subset(i2_table_per_replicate,i2_table_per_replicate$interval!="pos")
i3_table_per_replicate_forstat=subset(i3_table_per_replicate,i3_table_per_replicate$interval!="pos")

statTest_neg = c()
for(i in unique(neg_table_per_replicate_forstat$construct)){
  p = pairwise.wilcox.test(as.numeric(paste(neg_table_per_replicate_forstat$events)), neg_table_per_replicate_forstat$construct, p.adjust.method = "none") # test is here
  l = melt(p$p.value)[,3]
  print(l)
  names(l) = paste(melt(p$p.value)[,2], melt(p$p.value)[,1], sep="-")
  l = l[which(is.na(l)==FALSE)]
  t = multcompLetters(l,  compare='<', threshold = 0.05, Letters = c(letters, LETTERS, "."))
  tmp = data.frame(
    construct = names(t$Letters),
    letter = as.character(t$Letters),
    c = rep(i, length(t$Letters))
  )
  statTest_neg = rbind(statTest_neg, tmp)
}
statTest_neg$construct <- factor(statTest_neg$construct, c("mock", "empty", "sh1s",
                                                           "sh2s","sh3s","sh4s","sh4smutSS",
                                                           "sh5s","sh6s","sh6smutSS","sh7s","sh7smutSS","sh8s",
                                                           "sh9s","sh10s","sh10smutSS", "sh11s", "sh12s", "sh13s", "sh13smutSS"))
statTest_neg

pairwise.wilcox.test(neg_table_per_replicate_forstat$events, neg_table_per_replicate_forstat$construct, p.adjust.method = "none")


statTest_i1 = c()
for(i in unique(i1_table_per_replicate_forstat$construct)){
  p = pairwise.wilcox.test(as.numeric(paste(i1_table_per_replicate_forstat$events)), i1_table_per_replicate_forstat$construct, p.adjust.method = "none") # test is here
  l = melt(p$p.value)[,3]
  print(l)
  names(l) = paste(melt(p$p.value)[,2], melt(p$p.value)[,1], sep="-")
  l = l[which(is.na(l)==FALSE)]
  t = multcompLetters(l,  compare='<', threshold = 0.05, Letters = c(letters, LETTERS, "."))
  tmp = data.frame(
    construct = names(t$Letters),
    letter = as.character(t$Letters),
    c = rep(i, length(t$Letters))
  )
  statTest_i1 = rbind(statTest_i1, tmp)
}
statTest_i1$construct <- factor(statTest_i1$construct, c("sh1s",
                                                         "sh2s","sh3s","sh4s","sh4smutSS",
                                                         "sh5s","sh6s","sh6smutSS","sh7s","sh7smutSS","sh8s",
                                                         "sh9s","sh10s","sh10smutSS", "sh11s", "sh12s", "sh13s", "sh13smutSS"))
statTest_i1

pairwise.wilcox.test(i1_table_per_replicate_forstat$events, i1_table_per_replicate_forstat$construct, p.adjust.method = "none")


statTest_i2 = c()
for(i in unique(i2_table_per_replicate_forstat$construct)){
  p = pairwise.wilcox.test(as.numeric(paste(i2_table_per_replicate_forstat$events)), i2_table_per_replicate_forstat$construct, p.adjust.method = "none") # test is here
  l = melt(p$p.value)[,3]
  print(l)
  names(l) = paste(melt(p$p.value)[,2], melt(p$p.value)[,1], sep="-")
  l = l[which(is.na(l)==FALSE)]
  t = multcompLetters(l,  compare='<', threshold = 0.05, Letters = c(letters, LETTERS, "."))
  tmp = data.frame(
    construct = names(t$Letters),
    letter = as.character(t$Letters),
    c = rep(i, length(t$Letters))
  )
  statTest_i2 = rbind(statTest_i2, tmp)
}
statTest_i2$construct <- factor(statTest_i2$construct, c("sh1s",
                                                         "sh2s","sh3s","sh4s","sh4smutSS",
                                                         "sh5s","sh6s","sh6smutSS","sh7s","sh7smutSS","sh8s",
                                                         "sh9s","sh10s","sh10smutSS", "sh11s", "sh12s", "sh13s", "sh13smutSS"))
statTest_i2
pairwise.wilcox.test(i2_table_per_replicate_forstat$events, i2_table_per_replicate_forstat$construct, p.adjust.method = "none")


statTest_i3 = c()
for(i in unique(i3_table_per_replicate_forstat$construct)){
  p = pairwise.wilcox.test(as.numeric(paste(i3_table_per_replicate_forstat$events)), i3_table_per_replicate_forstat$construct, p.adjust.method = "none") # test is here
  l = melt(p$p.value)[,3]
  print(l)
  names(l) = paste(melt(p$p.value)[,2], melt(p$p.value)[,1], sep="-")
  l = l[which(is.na(l)==FALSE)]
  t = multcompLetters(l,  compare='<', threshold = 0.05, Letters = c(letters, LETTERS, "."))
  tmp = data.frame(
    construct = names(t$Letters),
    letter = as.character(t$Letters),
    c = rep(i, length(t$Letters))
  )
  statTest_i3 = rbind(statTest_i3, tmp)
}
statTest_i3$construct <- factor(statTest_i3$construct, c("sh1s",
                                                         "sh2s","sh3s","sh4s","sh4smutSS",
                                                         "sh5s","sh6s","sh6smutSS","sh7s","sh7smutSS","sh8s",
                                                         "sh9s","sh10s","sh10smutSS", "sh11s", "sh12s", "sh13s", "sh13smutSS"))
statTest_i3
pairwise.wilcox.test(i3_table_per_replicate_forstat$events, i3_table_per_replicate_forstat$construct, p.adjust.method = "none")


                ###########################################################################################
                #######2way ANOVA to see the impact of replicate and/or construct on the prop of cells ####
                #######                 (apply also to previous analyses)                              ####
                ###########################################################################################

library("ggpubr")

prop_intervals_per_replicate_forgraph <- events_per_intervals_per_replicate

#first test including only neg and i1,i2,i3 in this order
prop_intervals_per_replicate_forgraph_m_e=subset(prop_intervals_per_replicate_forgraph, prop_intervals_per_replicate_forgraph$construct!="mock"&prop_intervals_per_replicate_forgraph$construct!="empty"&prop_intervals_per_replicate_forgraph$interval!="pos")
prop_intervals_per_replicate_forgraph_m_e$interval <- factor(prop_intervals_per_replicate_forgraph_m_e$interval, c("neg","i1","i2","i3"))
prop_intervals_per_replicate_forgraph_m_e$construct <- factor(prop_intervals_per_replicate_forgraph_m_e$construct, c("sh1s",
                                                         "sh2s","sh3s","sh4s","sh4smutSS",
                                                         "sh5s","sh6s","sh6smutSS","sh7s","sh7smutSS","sh8s",
                                                         "sh9s","sh10s","sh10smutSS", "sh11s", "sh12s", "sh13s", "sh13smutSS"))

ggboxplot(prop_intervals_per_replicate_forgraph_m_e, x = "construct", y = "prop", color = "interval",fill="interval" ,palette = c("#666666","#1B9E77","#80c861","#66A61E"))
#ex:ggboxplot(my_data, x = "dose", y = "len", color = "supp", palette = c("#00AFBB", "#E7B800"))

#to do a propper anova, i need to remove alex's series and mock and empty
prop_intervals_per_replicate_forgraph_m_e = subset(prop_intervals_per_replicate_forgraph_m_e, prop_intervals_per_replicate_forgraph_m_e$replicate!="r6"&
                                                      prop_intervals_per_replicate_forgraph_m_e$replicate!="r7"&
                                                      prop_intervals_per_replicate_forgraph_m_e$replicate!="r8"&
                                                      prop_intervals_per_replicate_forgraph_m_e$replicate!="r9")

#test just en gardant neg et pos
prop_intervals_per_replicate_forgraph= subset(prop_intervals_per_replicate_forgraph, prop_intervals_per_replicate_forgraph$interval!="i1"
                                                           &prop_intervals_per_replicate_forgraph$interval!="i2"
                                                           &prop_intervals_per_replicate_forgraph$interval!="i3")


res.aov2 <- aov(events ~ replicate + construct + replicate:construct, data = prop_intervals_per_replicate_forgraph)
summary(res.aov2)
coefficients(res.aov2)
TukeyHSD(res.aov2, which = "construct")
#nothing comes out as significantly impacted by either replicate or constrcut

#no test w/o the neg interval
prop_intervals_per_replicate_forgraph_m_e=subset(prop_intervals_per_replicate_forgraph_m_e,prop_intervals_per_replicate_forgraph_m_e$interval!="neg")
res.aov2 <- aov(events ~ replicate + construct + replicate:construct, data = prop_intervals_per_replicate_forgraph_m_e)
summary(res.aov2)
coefficients(res.aov2)
TukeyHSD(res.aov2, which = "replicate")
TukeyHSD(res.aov2, which = "construct")

#still nothing comes out...


#test en rajoutant les replicat de alex mais en virant les mutSS
prop_intervals_per_replicate_forgraph=subset(prop_intervals_per_replicate_forgraph,prop_intervals_per_replicate_forgraph$interval!="pos")
prop_intervals_per_replicate_forgraph=subset(prop_intervals_per_replicate_forgraph,prop_intervals_per_replicate_forgraph$construct!="sh4smutSS"
                                             &prop_intervals_per_replicate_forgraph$construct!="sh6smutSS"
                                             &prop_intervals_per_replicate_forgraph$construct!="sh7smutSS"
                                             &prop_intervals_per_replicate_forgraph$construct!="sh10smutSS"
                                             &prop_intervals_per_replicate_forgraph$construct!="sh13smutSS"
                                             &prop_intervals_per_replicate_forgraph$construct!="mock")

res.aov2 <- aov(events ~ replicate + construct + replicate:construct, data = prop_intervals_per_replicate_forgraph)
summary(res.aov2)
coefficients(res.aov2)
TukeyHSD(res.aov2, which = "replicate")
TukeyHSD(res.aov2, which = "construct")

###############
# Contingency #
###############

library("FactoMineR")
library("factoextra")
library("gplots")
library(ggplot2)
library(ggpubr)

#BUILT TABLE
head(prop_intervals_per_replicate_forgraph)
table_to_transform_for_contingency=c()
for(x in unique(prop_intervals_per_replicate_forgraph$construct)){
  y = subset(prop_intervals_per_replicate_forgraph,prop_intervals_per_replicate_forgraph$construct==x)
  for (i in unique(y$interval)){
    j = subset(y,y$interval==i)
    events = sum(j$events)
    tmp = data.frame(construct=x, interval = i, events=events)
    tmp
    table_to_transform_for_contingency = rbind(table_to_transform_for_contingency, tmp)
  }
}
table_to_transform_for_contingency=subset(table_to_transform_for_contingency,table_to_transform_for_contingency$interval!="pos"&table_to_transform_for_contingency$construct!="mock"&table_to_transform_for_contingency$construct!="empty")
table_to_transform_for_contingency

#to get one column for each interval
library(tidyr)
dt<- pivot_wider(table_to_transform_for_contingency, names_from = interval, values_from = events) 
dt=data.frame(dt) # here is contingency table


#TEST KHI2
row.names(dt)=dt[,1]
dt
dt=dt[,-c(1)]
dt$neg=as.numeric(paste(dt$neg))
dt$i1=as.numeric(paste(dt$i1))
dt$i2=as.numeric(paste(dt$i2))
dt$i3=as.numeric(paste(dt$i3))
dt$i4=as.numeric(paste(dt$i4))
dt

khideux<-chisq.test(dt)
khideux
khideux$residuals

#BALLOONPLOT
#ggballoonplot(dt, fill = "value")+
# scale_fill_viridis_c(option = "C")+
# scale_color_viridis_c(option = "C")

###############
#PCA###########
###############
library(FactoMineR) #can't properly install the package


col=c("black","#666666","#1B9E77","#D95F02","#7570B3","#E7298A","#b22d71", "#66A61E", "#E6AB02", "#ae881b","#bb4732","#a14433","#80c861","#008080","#960018","#712834","#ffe4c4","#3f2a32","#6082b6","#626e7f")


#test with the 18 constructs of intersts
col_pca=c("black","#666666","#1B9E77","#D95F02","#7570B3","#E7298A","#b22d71", "#66A61E", "#E6AB02", "#ae881b","#bb4732","#a14433","#80c861","#008080","#960018","#712834","#ffe4c4","#3f2a32","#6082b6","#626e7f")
res.pca <- PCA(dt, graph = FALSE)
eig.val <- get_eigenvalue(res.pca)
eig.val
var <- get_pca_var(res.pca)
var
fviz_pca_biplot(res.pca, repel = TRUE,
                habillage ="none",
                geom.ind = "point",
                pointshape = 19,
                pointsize = 2.5,
                palette=col_pca,
                col.var = "black", 
                col.ind = c("Empty","1: Common", "2: Common - GC rich", "3: Common - AT rich","4: Rare", "5: Rare - GC rich", "6: Rare - AT rich")
                )

pca_fluo<-read.table("/media/marion/Commun/Documents/D2paper/Synthese_D2/R_codes_forGithub/Cytometry_Analyses/output/table_for_pca")
dt=pca_fluo[,-c(6)]
dt
res.pca <- PCA(dt)
eig.val <- get_eigenvalue(res.pca)
eig.val
var <- get_pca_var(res.pca)
var
fviz_pca_biplot(res.pca, repel = TRUE,
                habillage ="none",
                geom.ind = "point",
                pointshape = 19,
                pointsize = 3,
                palette=col_pca,
                col.var = "black", 
                col.ind = c("Empty","1: Common", "2: Common - GC rich", "3: Common - AT rich","4: Rare", "5: Rare - GC rich", "6: Rare - AT rich")
)
