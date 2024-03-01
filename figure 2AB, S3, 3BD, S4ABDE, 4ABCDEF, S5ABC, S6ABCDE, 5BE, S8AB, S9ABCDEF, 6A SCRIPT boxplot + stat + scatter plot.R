############################################################################################
#Script to realize composed graph with a linear regression plot and two adjacent boxplot   #
#showing the data spread of each variable.                                                 #
#original:marion picard,                                                                   #
#last update philippe paget 14/02/2024                                                     #
############################################################################################

#this script was used to generate figures 
#2AB, S3, 3BD, S4ABDE, 4ABCDEF, S5ABC, S6ABCDE, 5BE, S8AB, S9ABCDEF, 6A

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
library(cowplot)
####################
# Common knowledge
####################

#ALREADY DEFINE IN YOUR .CSV YOUR COLUMN NAMES THE WAY YOU WANT TO NAME YOUR AXIS 
#as such you can modify the entire script with "replace"

# avoid at all cost special characters in your dataframe, variable , line or column name 
# such as . , - * / \ any accent % $
# don"t include space either
# to seperate words in a file name or variable use _ (but if you can also avoid them it's better)

###
# a) import your df with push-buton system (avoid to ruin your header/ special characters during importation)

###
# b) change the name of your .csv file into what you want (will serve for the entire script)
big_recap <- beholdnorm2

###
# c) TRANSFORM YOUR columns containing experimental DATA INTO NUM OR INTEGER
big_recap[,10:83] <- lapply(big_recap[,10:83], as.numeric)
#if you add column to your .csv, don't forget to update the c) command accordingly


###
# d) you may need to select what part of your dataframe you want to analyze

#typically, we don't include mock and empty conditions in our graph
big_recap <- subset(big_recap, big_recap$construct!="mock")
big_recap <- subset(big_recap, big_recap$construct!="empty")

#you can also analyze sepcific replicate 
big_recap <- subset(big_recap, big_recap$replicate=="r1")

###
# e) choose the label and the order you want your conditions to appear on your graph, label and order (factor) have to match if you want correct graph plotting !, choose a color palette for your graphs
#ALL
big_recap$construct <- factor(big_recap$construct, c("empty","sh1","sh2","sh3","sh4sp","sh4","sh5","sh6sp","sh6","sh7sp","sh7","sh8","sh9","sh10sp","sh10", "sh11", "sh12", "sh13sp", "sh13"))
label= c("sh1","sh2","sh3","sh4sp","sh4","sh5","sh6sp","sh6","sh7sp","sh7","sh8","sh9","sh10sp","sh10", "sh11", "sh12", "sh13sp", "sh13")
col= c("#D95F02","#1B9E77","#6082b6","#E6AB02", "#ae881b", "#909495","#bb4732","#a14433","#80c861","#66A61E","#332298", "#CB2030","#E7298A","#b22d71","#3f2a32","#d1ab60","#6082b6","#7570B3")

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

#you can also download packages that do that automatically such as viridis


###
# g) Print any .csv you want
mypath="paste-brut-texte-here and replace anti-slash\\\\\ with slash ////, end with a slash"
#if you don't end it with a slash, it will save it one folder prior
write.table(df,file=paste(mypath,"XXX",".csv",sep=""),col.names = TRUE,row.names = FALSE)



                           #######################
                           # FIRST STAT TEST XXX #
                           #######################

#for clarity sake, think before going into R what you want your graph axis to be
# and already set that as column names in your .csv

# Replace XXX by the column name you want to analyze

statTest_X = c()
for(i in unique(big_recap$construct)){
  p = pairwise.wilcox.test(as.numeric(paste(big_recap$XXX)), big_recap$construct,paired=F, p.adjust.method = "BH") # test is here
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
  statTest_X = rbind(statTest_X, tmp)
}

#the A to Z filtering always put 10 before 1 so you need to reorder the construct to start w/ 1,2,3....viva l'alg?rie
statTest_X$construct <- factor(statTest_X$construct, label)
statTest_X


                ###################
                # The FIRST graph #
                ###################

#tiff("boxplot - XXX.tif", units="in", width=8, height=6, res=300, compression = 'none')

p_X <- ggplot(big_recap, aes(x=construct, y= XXX, fill = construct))+theme_ipsum()

p_X + geom_boxplot(width=0.8)+
  theme(legend.position="none")+
  geom_jitter(width = 0.15, size=5)+
  #ylim(0,10)+
  ylim(min(big_recap$XXX)-1,max(big_recap$XXX)+1)+
  #scale_y_continuous(trans = "log10")+     #transform les data en log10 mais uniquemment pour la rpz graphique
  scale_x_discrete(labels=label)+
  scale_fill_manual(values=col,name = "shble version", labels=label)+
  scale_color_manual(values=col,name = "shble version" , labels=label)+
  theme(axis.text.x = element_text(size=16, face ="bold", color = "black", angle=45),
        axis.text.y = element_text(size=16, face ="bold", color = "black"),
        axis.title.x = element_text(face="bold",size=18, hjust=0.5),
        axis.title.y = element_text(face="bold",size=18, hjust=0.5))+
  #ggtitle("") +
  #xlab("") + 
  ylab("XXX")
  geom_text(data = statTest_X, aes(x=construct, y=1+max(big_recap$XXX),label=statTest_X$letter),size=6) 
                                # in aes, add angle=90 tu turn the letters when too many

#dev.off()


                           #######################
                           #SECOND STAT TEST YYY #
                           #######################

# Replace YYY by the column name you want to analyze

statTest_Y = c()
for(i in unique(big_recap$construct)){
  p = pairwise.wilcox.test(as.numeric(paste(big_recap$YYY)), big_recap$construct, p.adjust.method = "BH") # test is here
  l = melt(p$p.value)[,3]
  print(l)
  names(l) = paste(melt(p$p.value)[,2], melt(p$p.value)[,1], sep="-")
  l = l[which(is.na(l)==F)]
  pvalues_Y_woNA <- data.frame(l)
  t = multcompLetters(l,  compare='<', threshold = 0.05, Letters = c(letters, LETTERS, "."))
  tmp = data.frame(
    construct = names(t$Letters),
    letter = as.character(t$Letters),
    c = rep(i, length(t$Letters))
  )
  statTest_Y = rbind(statTest_Y, tmp)
}

statTest_Y$construct <- factor(statTest_Y$construct, label)
statTest_Y


                    ####################
                    # The SECOND graph #
                    ####################

#tiff("boxplot - YYY.tif", units="in", width=6, height=4, res=300, compression = 'none')

p_Y <- ggplot(big_recap, aes(x=construct, y= YYY, fill=construct))+theme_ipsum()

p_Y + geom_boxplot(width=0.8)+
  theme(legend.position="none")+
  geom_jitter(width = 0.1, size=3)+
  #ylim(0,2.5)+
  #scale_y_continuous(trans = "log10")+     #transform les data en log10 mais uniquemment pour la rpz graphique
  scale_x_discrete(labels=label)+
  scale_fill_manual(values=col,name = "shble version", labels=label)+
  scale_color_manual(values=col,name = "shble version" , labels=label)+
  theme(plot.title= element_text(face="bold",size=16), axis.text.x = element_text(size=12, color = "black"),axis.text.y = element_text(size=14, color = "black"),axis.title.x = element_text(face="bold",size=14, hjust=0.5),axis.title.y = element_text(face="bold",size=12, hjust=0.5))+
  ggtitle("") +
  xlab("") + 
  ylab("YYY")+
  geom_text(data = statTest_Y, aes(x=construct, y=max(big_recap$XXX)+1, label=statTest_Y$letter),size=6) 
                                # in aes, add angle=90 tu turn the letters when too many

#dev.off()

                                #############
                                ##  X vs Y ##
                                #############

#be very cautious with what is displayed on graph
#because adding the stat_cor sometimes create a mismatch between color and data !

###################################################################################
#always generate a version of your graph when you segregate the data by replicate 
# to spot potential biais between experiments
#label= c("r1","r2","r3","r4","r5","r6","r7","r8")
#big_recap$replicate<- factor(big_recap$replicate, c("r1","r2","r3","r4","r5","r6","r7","r8"))
#col= c("#D95F02","#1B9E77","#6082b6","#E6AB02","#909495","#bb4732","#80c861","#332298","#CB2030","#3f2a32","#7570B3")
###################################################################################

###################################################
###careful with the limit of the graph cause if datapoints are outside of graphs boundaries they are not taken into account for lm, pvalue and R coeff calculation 

#tiff("scatter - XXX vs. YYY.tif", units="in", width=6, height=4, res=300, compression = 'none')

sp <-ggscatter(big_recap, x ="XXX" , y ="YYY" ,na.rm=T, 
                add = "reg.line", #correlation of each condition defined by color=""   #other fit option = "loess" instead of "reg.line"
                conf.int =F,
                size = 5, alpha = 0.6,
                color="construct")+
 # geom_smooth(method=lm, fill = "gray", fullrange=T, color = "black")+ #formula=y~x+0 #general correlation of all the data
  scale_fill_manual(values=col,name = "version", labels = label)+
  scale_color_manual(values=col,name = "version", labels = label)+ 
  theme(legend.position = "bottom")+
  xlab("XXX") + 
  ylab("YYY")+  
  theme(axis.text.x = element_text(size=16, face ="bold", color = "black"),
        axis.text.y = element_text(size=16, face ="bold", color = "black"),
        axis.title.x = element_text(face="bold",size=18, hjust=0.5),
        axis.title.y = element_text(face="bold",size=18, hjust=0.5))+
  ggtitle("") +
  guides(color = guide_legend(nrow = 2)) #on how many line do you want your legend to appear

#here you add the correlation p value and coeff    
#of the general correlation
sp2 <-sp + stat_cor(label.y.npc = 0.9,label.x.npc = 0.5,size=6, method="pearson" ) #label.x.npc = "right" #aes(color= cub),
sp2

#of each condition defined by color=
sp2 <-sp + stat_cor(aes(color=construct),label.y.npc = 1,label.x.npc = 0.6,size=6, method="pearson" ) #label.x.npc = "right"
sp2

#to display regression line equation
spslope <- sp + stat_regline_equation(aes(color= construct),label.y.npc = 1,label.x.npc = 0.7)
spslope

#to add boxplot for x and Y on either side ofthe graph
ggMarginal(sp2, type = "boxplot")

#dev.off()

##test if reg line intersect is stat different from the origin
lm <- lm(formula = ibaq_shble ~ ibaq_gfp, data= big_recap)
summary(lm)
lm0 <- lm(formula = ibaq_shble ~ 0 + ibaq_gfp, data= big_recap)
summary(lm0)
anova(lm, lm0)

