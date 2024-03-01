#####################################################
#            TEST PCA data shbleSTOP                #
#####################################################

library(janitor)
library(ade4)
library(factoextra)
library(FactoMineR)
library(magrittr)
library(ggplot2)
library(cowplot)
library(ggpubr)

#put the first row as header (moins chiant que read .csv)
cousin <- input_PCA_arn_ibaq
cousin <- subset(cousin, cousin$cub!="no")
cousin <- cousin[,-2:-6]
rownames(cousin) <- cousin[,1]
cousin <- cousin[,-1]
cousin <- cousin[,-1:-4]
#converti en valeur num?rique
cousin[,1:3] <- lapply(cousin[,1:3], as.numeric)


##########
#la PCA###
##########

####en utilisant le folding de 5'utr+CDSshble
cousin_foldingCDS.pca <- dudi.pca(cousin,
                           scannf = T,  
                           nf = 5
)


plot_weigthdimension <- fviz_eig(cousin_foldingCDS.pca)

plot_PCA <-fviz_pca_ind(cousin_foldingCDS.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
plot_PCA

plot_vecteurvariable <-fviz_pca_var(cousin_foldingCDS.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

biplotpca <-fviz_pca_biplot(cousin_foldingCDS.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
biplotpca

plot_grid(plot_weigthdimension, plot_PCA, plot_vecteurvariable, biplotpca, labels = "AUTO")

################r?cup?ration des coordonn?es des axes de la pca
# aller dans df cousin.pca
  #li
    # df axis 1, axis 2, ...

coor_dimension_pcaCDS <- cousin_foldingCDS.pca[["li"]]

#cousin_nomconstruct <- cousin[,1]

coor_dimension_pcaCDS <- cbind(cousin, coor_dimension_pcaCDS)

# Print big_recap
mypath="XXX"
write.table(coor_dimension_pcaCDS,file=paste(mypath,"dimension_pca",".csv",sep=""),col.names = TRUE,row.names = FALSE)
