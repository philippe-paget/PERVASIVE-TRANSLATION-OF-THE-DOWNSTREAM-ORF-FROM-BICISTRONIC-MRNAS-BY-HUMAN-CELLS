############################################################################################
#Script to realize composed graph with a linear regression plot and two adjacent boxplot   #
#showing the data spread of each variable.                                                 #
#original:marion picard,                                                                   #
#last update philippe paget 14/02/2024                                                     #
############################################################################################

#this script was used to generate figures 
#3A, 3C, 5C and 5D

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

# avoid at all cost special characters in your dataframe, variable , line or column name 
# such as . , - * / \ any accent % $
# don"t include space either
# to seperate words in a file name or variable use _ (but if you can also avoid them it's better)

###
# a) import your df with push-buton system (avoid to ruin your header/ special characters during importation)

###
# b) change the name of your .csv file into what you want (will serve for the entire script)
big_recap <- boxplot_diffSHBLEtoGFPiBAQratio

###
# c) TRANSFORM YOUR columns containing experimental DATA INTO NUM OR INTEGER

big_recap[,9:11] <- lapply(big_recap[,9:11], as.numeric)
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

###########################################################################

big_recap$construct <- factor(big_recap$protein, c("SHBLE", "GFP", "GFP_mono"))
label = c("SHBLE", "GFP", "GFP_mono")

#selection for pairwise comparison
#GFP vs.GFP_mono
big_recap <- subset(big_recap, big_recap$protein!="SHBLE")
montest=wilcox.test(value~protein, data=big_recap, exact=T)
montest

#GFP vs. SHBLE
big_recap <- subset(big_recap, big_recap$protein!="GFP_mono")
montest=wilcox.test(value~protein, data=big_recap, exact=F)
montest

#SHBLE vs. GFP_mono
big_recap <- subset(big_recap, big_recap$protein!="GFP")
montest=wilcox.test(value~protein, data=big_recap, exact=T)
montest

####################################################

big_recap_match <- subset(big_recap, big_recap$match=="match")
big_recap_under <- subset(big_recap, big_recap$match=="undermatch")
big_recap_over <- subset(big_recap, big_recap$match=="overmatch")

shapiro.test(big_recap_match$value)
shapiro.test(big_recap_under$value)
shapiro.test(big_recap_over$value)
# all three are not normal so non-para wilcox test w/ continuity

statTest_X = c()
for(i in unique(big_recap$match)){
  p = pairwise.wilcox.test(as.numeric(paste(big_recap$value)), big_recap$match,paired=F, p.adjust.method = "BH") # test is here
  l = melt(p$p.value)[,3]
  print(l)
  names(l) = paste(melt(p$p.value)[,2], melt(p$p.value)[,1], sep="-")
  l = l[which(is.na(l)==F)]
  pvalues_stat_test_X_woNA <- data.frame(l) #added this to retreive the pvalues w/o NaN
  t = multcompLetters(l,  compare='<', threshold = 0.05, Letters = c(letters, LETTERS, "."))
  tmp = data.frame(
    match = names(t$Letters),
    letter = as.character(t$Letters),
    c = rep(i, length(t$Letters))
  )
  statTest_X = rbind(statTest_X, tmp)
}

#########################DOT-PLOT#################################
big_recap <- subset(big_recap, big_recap$protein!="ratio")
big_recap$match<- factor(big_recap$match, c("overmatch","match","undermatch"))
label <- c("sh1","sh2","sh3","sh4","sh5","sh6","sh7","sh8","sh9","sh10","sh11","sh12","sh13")
label2 <-c("overmatch","match","undermatch")
col= c("#D95F02","#1B9E77","#6082b6")

#tiff("SP - log10 SHBLEoverGFP levels mutated vs. spliced boxplot.tif", units="in", width=8, height=6, res=300, compression = 'none')

p_X <- ggplot(big_recap, aes(x=match, y= value , fill= construct))+theme_ipsum()

p_X + geom_boxplot(width=0.8)+
  theme(legend.position="bottom")+
  #ylim(0,2000)+
  scale_y_continuous(trans = "log10")+     #transform les data en log10 mais uniquemment pour la rpz graphique
  scale_x_discrete(labels=label2)+
  scale_fill_manual(values=col,name = "shble version", labels=label)+
  scale_color_manual(values=col,name = "shble version" , labels=label2)+
  xlab("") + 
  ylab("LOG10 SHBLE/GFP levels")+
  theme(axis.text.x = element_text(size=16, face ="bold", color = "black"),
        axis.text.y = element_text(size=16, face ="bold", color = "black"),
        axis.title.x = element_text(face="bold",size=18, hjust=0.5),
        axis.title.y = element_text(face="bold",size=18, hjust=0.5))
geom_text(data = statTest_X, aes(x=mut, y=6000,label=statTest_X$letter),size=6) # in aes, add angle=90 tu turn the letters when too many
geom_text(data = pvalues_stat_test_X_woNA, aes(x=match, y=6000,label=statTest_X$letter),size=6)

#dev.off()

p_X + geom_dotplot(binaxis = "y",stackdir='center',  binwidth = 0.2,
                   position=position_dodge(1))+
  theme(legend.position="bottom")+
  #ylim(0,4000000000000)+
  scale_y_continuous(trans = "log10")+     #transform les data en log10 mais uniquemment pour la rpz graphique
  scale_x_discrete(labels=label)+
  scale_fill_manual(values=col,name = "protein", labels=label2)+
  scale_color_manual(values=col,name = "protein" , labels=label2)+
  xlab("") + 
  ylab("LOG10 Protein levels (norm. iBAQ)")+
  theme(axis.text.x = element_text(size=16, face ="bold", color = "black"),
        axis.text.y = element_text(size=16, face ="bold", color = "black"),
        axis.title.x = element_text(face="bold",size=18, hjust=0.5),
        axis.title.y = element_text(face="bold",size=18, hjust=0.5))

#dev.off()

##############################################################################

#t######test ratio splicing
big_recap <- subset(big_recap, big_recap$splicing=="yes")
#selec for comparaison deux à deux
big_recap <- subset(big_recap, big_recap$protein=="ratio")
montest=wilcox.test(value~mut, data=big_recap, exact=F)
montest