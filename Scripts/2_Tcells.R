############################################
############################################
## Analysis of efficaccy - Takeda Vaccine ##
############################################
############### FLOW DATA ##################
################ T Cells ###################
############################################

library(readxl)
library(PCAtools)
library(factoextra)
library(FactoMineR)
library(Rtsne)
library(ggrepel)
library(clustree)
library(ComplexHeatmap)
library(pals)
library(RColorBrewer)
library(viridis)
library(fmsb)


set.seed(123)

tcells <- read_xlsx("Data/DENV_T_cell_summary_CLEAN 2.xlsx")
baseline <- tcells$`ELISA V1 positive?`
for(i in 1:length(baseline)){
  if(is.na(baseline[i])){
    baseline[i] <- "Negative"
  }
}
sex <- tcells$Sex
Age <- tcells$Age
for(i in 1:length(Age)){
  if(Age[i] < 65){
    Age[i] <- 1
  } else if(Age[i] >= 65){
    Age[i] <- 2
  }
}

visist <- tcells$`Time points`



# Splitting into the timepoints

tcells_V1 <- subset(tcells, subset = `Time points` == "V1")

tcells_V3 <- subset(tcells, subset = `Time points` == "V3")

#### PCA total

tcells_bio <- tcells[,9:42]#38


pca_total <- prcomp(scale(na.omit(tcells_bio)))

fviz_pca_ind(pca_total, label = "none")


# How many PCs should we use?

fviz_screeplot(pca_total, addlabels = TRUE, ylim = c(0, 50))

fviz_pca_var(pca_total, col.var = "black")

# There are clearly two types of celular response

pdf("Figures/PCA_tcells.pdf")

fviz_pca_ind(pca_total, label = "none", addEllipses = T)

dev.off()


# Is id due to the previous infection?

pdf("Figures/PCA_tcells_baseline.pdf")

fviz_pca_ind(pca_total, habillage = baseline,  label = "none", palette = c("#18AFD3", "#840EC9"), addEllipses = T)

dev.off()


# Is id due to the time after injection?

pdf("Figures/PCA_tcells_visits.pdf")

fviz_pca_ind(pca_total, habillage = visist,  label = "none", palette = c("#B15A28", "#4DAF4A"), addEllipses = T)

dev.off()


# Is id due to the sex?

pdf("Figures/PCA_tcells_Sex.pdf")

fviz_pca_ind(pca_total, habillage = sex,  label = "none", palette = c("#67ADDF", "#102B43"), addEllipses = T)

dev.off()


# Is id due to Age?

pdf("Figures/PCA_tcells_Age.pdf")

fviz_pca_ind(pca_total, habillage = Age,  label = "none", palette = c("#67ADDF", "#102B43"), addEllipses = T)

dev.off()



#### PCA specific


tcells_V1 <- tcells_V1[,9:42]

pca <- prcomp(na.omit(scale(tcells_V1)))


pdf("Figures/PCA_tcells_V1.pdf")
fviz_pca_ind(pca)
dev.off()


tcells_V3 <- tcells_V3[,9:42]

pca <- prcomp(na.omit(scale(tcells_V3)))

pdf("Figures/PCA_tcells_V3.pdf")
fviz_pca_ind(pca)
dev.off()



# K means clustering

km <- kmeans(as.matrix(scale(tcells_bio)), centers = 2)

table(km[[1]])


pdf("Figures/kmeans_clustering.pdf")
fviz_pca_ind(pca_total, habillage = km[[1]])
dev.off()


