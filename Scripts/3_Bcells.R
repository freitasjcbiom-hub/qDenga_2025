###########################################
###########################################
## Analysis of eficaccy - Takeda Vaccine ##
###########################################
############### FLOW DATA #################
################ B Cells ##################
###########################################


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
library(NLP)
library(randomForest)
library(janitor)
library(Hmisc)
library(corrplot)
library(viridis) 


# Setting a global seed

set.seed(123)



#Starting the code



bcells <- read_xlsx("Data/B_cell_data_summary_LC_CLEAN.xlsx", sheet = 3)

colnames(bcells)

pca_total <- prcomp(na.omit(bcells[,3:ncol(bcells)]))

pdf("Figures/CleanBCells.pdf")

fviz_pca_ind(pca_total)

dev.off()



pdf("Figures/Bcells_Timepoints.pdf")

fviz_pca_ind(pca_total, habillage = na.omit(bcells)$Timepoints)

dev.off()
 #PCA shows a transition (from right to left) of the first timepoint to the last one



# V1 x Non-V1

new_data <-  na.omit(bcells)

Timepoints_2 <- c()

for(i in 1:nrow(new_data)){
  if(new_data$Timepoints[i] == "V1"){
    Timepoints_2[i] <- "V1"
  } else if(new_data$Timepoints[i] != "V1"){
    Timepoints_2[i] <- "Non-V1"
  }
}


pdf("Figures/Bcells_Timepoints_2.pdf")

fviz_pca_ind(pca_total, habillage = Timepoints_2)

dev.off()


# TSNE

data_tsne <- na.omit(bcells[,3:ncol(bcells)])


set.seed(6)
Rtsne <- Rtsne(data_tsne, perplexity=5, check_duplicates = FALSE)
rtsne_df <- as.data.frame(Rtsne$Y)
rtsne_df$names <- new_data$`Patient ID`


saveRDS(rtsne_df, "Figures/bcell_tsne.rds")

rtsne_df <- readRDS("/Users/jramalh2/Documents/Yale/bcell_tsne.rds")


rtsne_plot <- ggplot(rtsne_df, aes(x=V1, y=V2)) + 
  geom_point( colour = "black", shape = 21, size = 4, stroke = 1) + 
  scale_fill_brewer(palette = "Set2") + 
  geom_text_repel(aes(label=rownames(rtsne_df)), size=4) +
  theme_minimal()



rtsne <- ggplot(rtsne_df, aes(x=V1, y=V2)) + 
  geom_point(aes(fill=new_data$Timepoints), colour = "black", shape = 21, size = 4, stroke = 1) + 
  theme_classic()

pdf("Figures/tSNE_Timepoints.pdf")

rtsne

dev.off()




km <- kmeans(data_tsne, centers = 2)




rtsne_HC5 <- ggplot(rtsne_df, aes(x=V1, y=V2)) + 
  geom_point(aes(fill=as.character(km[[1]])), colour = "black", shape = 21, size = 4, stroke = 1) + 
  theme_classic()


pdf("Figures/tSNE_kMeans.pdf")

rtsne_HC5

dev.off()



clusters <- data.frame(km[[1]],  na.omit(bcells)$Timepoints)

colnames(clusters) <- c("Km", "Timepoint")

clusters_1 <- clusters[clusters$Km == 1, ]
clusters_2 <- clusters[clusters$Km == 2, ]


table(clusters_1$Timepoint)
table(clusters_2$Timepoint)


#### Identifying patients already in cluster 1 since V1


bcells_2 <- na.omit(bcells)
bcells_2$cluster <- km[[1]]

# Crossing with the demographic information

tsne_demographics_3 <- readRDS("RDS_objects/tsne_demographics_3.rds")

bcells_2$ID <- bcells_2$`Patient ID`

bcells_2_demo <- merge(bcells_2, tsne_demographics_3, by = 'ID')



# Subsetting on V1 data:

V1_bcells <- bcells_2_demo[bcells_2_demo$Timepoints == "V1", ]

table(V1_bcells$cluster)


V1_bcells_cluster1 <- V1_bcells[V1_bcells$cluster == 1, ]

table(V1_bcells_cluster1$hclust_5)








table(clusters_1$re)


