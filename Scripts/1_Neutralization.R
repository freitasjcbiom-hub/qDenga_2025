###########################################
###########################################
## Analysis of eficaccy - Takeda Vaccine ##
###########################################
########## NEUTRALIZATION DATA ############
###########################################
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
library(randomForest)
library(treemap)
library(tidyverse)


###### Setting general seed
set.seed(123)



####------------ Data Loading

data <- read_xlsx("Data/Meta_table.xlsx")


####------------ Data Cleaning


source("Functions/cleanData.R")


# Excluding the rows with more than 10 NA and the collumns with more than 20 NA
data_filtered <-  cleanData(data, 20, 10)[[1]]



# Organizing the data


ID <- data_filtered$ID

ID_2 <- paste0("Patient_", ID)


source("Functions/df_to_numMatrix.R")

data_filtered_num <- as.data.frame(df_to_numMatrix(data_filtered))

rownames(data_filtered_num) <- ID_2



# Creating a collum that splits the data between young and old people

age_2 <- c()
for(i in 1:nrow(data_filtered_num)){
  if(data_filtered_num$Age[i] > 50){
    age_2 <- c(age_2, 2)
  } else if(data_filtered_num$Age[i] < 50){
    age_2 <- c(age_2, 1)
  } else if(is.na(data_filtered_num$Age[i])){
    age_2 <- c(age_2, '')
  }
}

data_filtered_num$age_2 <- as.character(age_2)


### Splitting data between demographic and response data

demographic <- cbind(data_filtered_num[ ,1:8], Age_2 = data_filtered_num[ ,25])
demographic$IDs <- rownames(demographic)

vac_res <- data_filtered_num[ ,9:24]






#####------------ Data analysis





##-----------PCA


# Patients PCA


pca <- prcomp(na.omit(vac_res))

fviz_pca_ind(pca)


pdf("Figures/First_PCA.pdf")

fviz_pca_ind(pca)

dev.off()



##---------- Hirarchical clustering


dd <-  dist(vac_res)

hc <- hclust(dd)
plot(hc)


pdf("Figures/Hierarchical_clustering.pdf")
plot(hc)
dev.off()



cluster <- as.data.frame(cutree(hc, k=4))
cluster$IDs <- rownames(cluster)



demographic <- merge(demographic, cluster, by = "IDs")
rownames(demographic) <- demographic$IDs


### Tsne

set.seed(4231)
Rtsne <- Rtsne(na.omit(vac_res), perplexity=5, check_duplicates = FALSE)
rtsne_df <- as.data.frame(Rtsne$Y)
rtsne_df$names <- rownames(na.omit(vac_res))


rtsne_plot <- ggplot(rtsne_df, aes(x=V1, y=V2)) + 
  geom_point( colour = "black", shape = 21, size = 4, stroke = 1) + 
  scale_fill_brewer(palette = "Set2") + 
  geom_text_repel(aes(label=rownames(rtsne_df)), size=4) +
  theme_minimal()




# Clustering

hc.norm = hclust(dist(Rtsne$Y))
rtsne_df$hclust_1 = factor(cutree(hc.norm, 1))
rtsne_df$hclust_2 = factor(cutree(hc.norm, 2))
rtsne_df$hclust_3 = factor(cutree(hc.norm, 3))
rtsne_df$hclust_4 = factor(cutree(hc.norm, 4))
rtsne_df$hclust_5 = factor(cutree(hc.norm, 5))







#saveRDS(rtsne_df, "RDS_objects//tsne.rds")


rtsne_HC5 <- ggplot(rtsne_df, aes(x=V1, y=V2)) + 
  geom_point(aes(fill=hclust_5), colour = "black", shape = 21, size = 4, stroke = 1) + 
  theme_classic()
rtsne_HC5


pdf("Figures/tSNE_2.pdf")

plot(rtsne_HC5)

dev.off()






# Adding the demographic data to the Tsne matrix

demographic$names <- rownames(demographic)

tsne_demographics <- merge(rtsne_df, demographic, by="names")


saveRDS(tsne_demographics, "RDS_objects/tsne_demographics.rds")

#tsne_demographics <- readRDS("RDS_objects/tsne_demographics.rds")

vac_res$names <- rownames(vac_res)

tsne_demographics_2 <- merge(tsne_demographics, vac_res, by = "names")


write.csv(tsne_demographics_2, "Data/demographics.csv")


# Checking the pattern of serostatus baseline


tsne_demographics$`Serostatus Baseline` <- as.factor(tsne_demographics$`Serostatus Baseline`)
rtsne_SeroBase <- ggplot(tsne_demographics, aes(x=V1, y=V2)) + 
  geom_point(aes(fill=`Serostatus Baseline`), colour = "black", shape = 21, size = 4, stroke = 1) + 
  theme_classic()


pdf("Figures/tSNE_serostatus.pdf")
plot(rtsne_SeroBase)
dev.off()

# Checking the pattern of Yellow Fever Vaccine 

tsne_demographics$`Vacina Febre Amarela` <- as.factor(tsne_demographics$`Vacina Febre Amarela`)
rtsne_FA <- ggplot(tsne_demographics, aes(x=V1, y=V2)) + 
  geom_point(aes(fill=`Vacina Febre Amarela`), colour = "black", shape = 21, size = 4, stroke = 1) + 
  theme_classic()


pdf("Figures/tSNE_FebreAmarela.pdf")
plot(rtsne_FA)
dev.off()


# Checking the pattern of Age

tsne_demographics$Age_2 <- as.factor(tsne_demographics$Age_2)
rtsne_Age <- ggplot(tsne_demographics, aes(x=V1, y=V2)) + 
  geom_point(aes(fill=`Age_2`), colour = "black", shape = 21, size = 4, stroke = 1) + 
  theme_classic()


pdf("Figures/tSNE_Age.pdf")
plot(rtsne_Age)
dev.off()



# Checking the pattern of Sex



tsne_demographics$Sex <- as.factor(tsne_demographics$Sex)
rtsne_Sex <- ggplot(tsne_demographics, aes(x=V1, y=V2)) + 
  geom_point(aes(fill=Sex), colour = "black", shape = 21, size = 4, stroke = 1) + 
  theme_classic()

pdf("Figures/tSNE_Sex.pdf")
plot(rtsne_Sex)
dev.off()



# Barplot of demographics through the clusters 


tsne_demographics$`Sex` <- as.character(tsne_demographics$`Sex`)
barplot <- ggplot(data = tsne_demographics) + 
  geom_bar(mapping = aes(x = hclust_5, fill =`Sex`), position = "fill") +
  theme_classic()

pdf("Figures/Barplot_Sex.pdf")
barplot
dev.off()


tsne_demographics$`Serostatus Baseline` <- as.character(tsne_demographics$`Serostatus Baseline`)
barplot_2 <- ggplot(data = tsne_demographics) + 
  geom_bar(mapping = aes(x = hclust_5, fill =`Serostatus Baseline`), position = "fill") +
  theme_classic()

pdf("Figures//Barplot_Baseline.pdf")
barplot_2
dev.off()


tsne_demographics$`Vacina Febre Amarela` <- as.character(tsne_demographics$`Vacina Febre Amarela`)
barplot_3 <- ggplot(data = tsne_demographics) + 
  geom_bar(mapping = aes(x = hclust_5, fill =`Vacina Febre Amarela`), position = "fill") +
  theme_classic()


pdf("Figures/Barplot_VFA.pdf")
barplot_3
dev.off()


tsne_demographics$Age_2 <- as.character(tsne_demographics$Age_2)
barplot_4 <- ggplot(data = tsne_demographics) + 
  geom_bar(mapping = aes(x = hclust_5, fill = Age_2), position = "fill") +
  theme_classic()


pdf("Figures/Barplot_Age2.pdf")
barplot_4
dev.off()



#--------------------------------- Heatmap



#taking clustering from tsne

clustering <- cbind(names = rtsne_df$names, hc_5 = rtsne_df$hclust_5)

#Merging it to neutralization data 

vac_res$names <- rownames(vac_res)

#Adjusting the data to the Heatmap

heatmap_data <- merge(clustering, vac_res, by = "names")

heatmap_data$names <- NULL

hc_5 <- heatmap_data$hc_5
heatmap_data$hc_5 <- NULL

scaled_heatmap_data <- scale(heatmap_data)

#Definying the rows split

denv_sorotipos <- c(rep("1", 4), rep("2", 4), rep("3", 4), rep("4", 4))

#Heatmap

HM <- ComplexHeatmap::Heatmap(t((scaled_heatmap_data)),
                              column_split = tsne_demographics$hclust_5, 
                              row_split = denv_sorotipos,
                              cluster_columns = T, 
                              cluster_rows = F, 
                              col= c("#abd9e9", "#74add1", "#4575b4","#8e44ad", "#840EC9"))


# Adapting the range to take care of possible ooutliers

source("Functions/newRange.R")

scaled_heatmap_data_InRange <- newRange(scaled_heatmap_data, 2.5, -2.5)

HM <- ComplexHeatmap::Heatmap(t((scaled_heatmap_data_InRange)),
                              column_split = tsne_demographics$hclust_5, 
                              row_split = denv_sorotipos,
                              cluster_columns = T, 
                              cluster_rows = F, 
                              col= c("#abd9e9", "#74add1", "#4575b4","#8e44ad", "#840EC9"))

draw(HM)


pdf("Figures//Heatmap_hc5_3.pdf")
plot(HM)
dev.off()


#--------------
# ISSUE: the normalization above is masking low personal increases in the neutralization 
# to adress it, we will perform another heatmap but with the fold change of each patient 
#-------------


#### Calculating the FC


# Taking out of the log


antilog_vac_res <- 10^as.data.frame(df_to_numMatrix(vac_res))

View(antilog_vac_res)


# FC

FC_table_1 <- data.frame(row.names = rownames(antilog_vac_res))

for(i in 1:4){
  FC_table_1 <- cbind(FC_table_1, antilog_vac_res[i]/antilog_vac_res[1])
}


FC_table_2 <- data.frame(row.names = rownames(antilog_vac_res))

for(i in 5:8){
  FC_table_2 <- cbind(FC_table_2, antilog_vac_res[i]/antilog_vac_res[5])
}


FC_table_3 <- data.frame(row.names = rownames(antilog_vac_res))

for(i in 9:12){
  FC_table_3 <- cbind(FC_table_3, antilog_vac_res[i]/antilog_vac_res[9])
}


FC_table_4 <- data.frame(row.names = rownames(antilog_vac_res))

for(i in 13:16){
  FC_table_4 <- cbind(FC_table_4, antilog_vac_res[i]/antilog_vac_res[13])
}


FC_df <- cbind(FC_table_1, c(FC_table_2, FC_table_3, FC_table_4))

col_names <- colnames(FC_df)

for(i in 1:length(col_names)){
  col_names[i] <- gsub('Log','',col_names[i])
}

colnames(FC_df) <- col_names

write.csv(FC_df, "Data/FC_Table.csv")

#FC_df <- read.csv("Data/FC_Table.csv")


# Creating a Heatmap with the FC


HM_FC_pos <- ComplexHeatmap::Heatmap(t(na.omit(FC_df)),
                                     column_split = tsne_demographics$hclust_5,
                                     row_split = denv_sorotipos,
                                     show_column_names = F,
                                     cluster_columns = T, 
                                     cluster_rows = F, 
                                     col = magma(10))


# Adjusting the range so we can se the heatmap in better detail

FC_df_2 <- na.omit(FC_df)

FC_df_2 <- newRange(FC_df_2, 10, min(FC_df_2))

HM_FC_pos <- ComplexHeatmap::Heatmap(t(FC_df_2),
                                     column_split = tsne_demographics$hclust_5,
                                     row_split = denv_sorotipos,
                                     show_column_names = F,
                                     cluster_columns = T, 
                                     cluster_rows = F, 
                                     col = magma(10))




HM_FC_pos <- draw(HM_FC_pos)

pdf("Figures/FC_Heatmap_MAX10.pdf")
HM_FC_pos
dev.off()


#--------------------------------- Breadth and Imprinting


### How many patients had a fold change higher than 3?

FC_df_2 <- as.data.frame(FC_df_2)

a_1 <- FC_df_2$`D1 - V2 FRNT50`>3
table(a_1)
sum(a_1)/nrow(FC_df_2)
pdf("Figures/pie_V2D1.pdf")
pie(table(a_1), col = brewer.set1(4), labels = c("Under Treshold", "Above Treshold"))
dev.off()


a_2 <- FC_df_2$`D2 - V2 FRNT50`>3
table(a_2)
sum(a_2)/nrow(FC_df_2)
pdf("Figures/pie_V2D2.pdf")
pie(table(a_2), col = brewer.set1(4), labels = c("Under Treshold", "Above Treshold"))
dev.off()



a_3 <- FC_df_2$`D3 - V2 FRNT50`>3
table(a_3)
sum(a_3)/nrow(FC_df_2)
pdf("Figures/pie_V2D3.pdf")
pie(table(a_3), col = brewer.set1(4), labels = c("Under Treshold", "Above Treshold"))
dev.off()




a_4 <- FC_df_2$`D4 - V2 FRNT50`>3
table(a_4)
sum(a_4)/nrow(FC_df_2)
pdf("Figures/pie_V2D4.pdf")
pie(table(a_4), col = brewer.set1(4), labels = c("Under Treshold", "Above Treshold"))
dev.off()




# Stacked barplots showing the distribution of patients above and bellow the treshold

tresholds_V2 <- list(a_1, a_2, a_3, a_4)


for(i in 1:length(tresholds_V2)){
  
  stack_bars <- as.data.frame(as.vector(tresholds_V2[[i]]))
  
  
  stack_bars$type <- rep(c(paste0("D", i)), nrow(stack_bars))
  colnames(stack_bars) <- c("Treshold", "D_type")
  
  barplot <- ggplot(data = stack_bars) + 
    geom_bar(mapping = aes(x = D_type, fill = Treshold), position = "fill") +
    scale_fill_manual("Legend", values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"))+
    theme_classic()
  
  assign("path_exp", paste0("Figures/Barplot_V2D", i, ".pdf"))
  
  pdf(path_exp)
  plot(barplot)
  dev.off()
  
}




# What about the thresholds in V4?



b_1 <-  FC_df_2$`D1 - V4 FRNT50`>3
table(b_1)
sum(b_1)/nrow(FC_df_2)
pdf("Figures/pie_V4D1.pdf")
pie(table(b_1), col = brewer.set1(4), labels = c("Under Treshold", "Above Treshold"))
dev.off()


b_2 <- FC_df_2$`D2 - V4 FRNT50`>3
table(b_2)
sum(b_2)/nrow(FC_df_2)
pdf("Figures/pie_V4D2.pdf")
pie(table(b_2), col = brewer.set1(4), labels = c("Under Treshold", "Above Treshold"))
dev.off()


b_3 <- FC_df_2$`D3 - V4 FRNT50`>3
table(b_3)
sum(b_3)/nrow(FC_df_2)
pdf("Figures/pie_V4D3.pdf")
pie(table(b_3), col = brewer.set1(4), labels = c("Under Treshold", "Above Treshold"))
dev.off()


b_4 <- FC_df_2$`D4 - V4 FRNT50`>3
table(b_4)
sum(b_4)/nrow(FC_df)
pdf("Figures/pie_V4D4.pdf")
pie(table(b_4), col = brewer.set1(4), labels = c("Under Treshold", "Above Treshold"))
dev.off()



tresholds_V4 <- list(b_1, b_2, b_3, b_4)


for(i in 1:length(tresholds_V4)){
  
  stack_bars <- as.data.frame(as.vector(tresholds_V4[[i]]))
  
  
  stack_bars$type <- rep(c(paste0("D", i)), nrow(stack_bars))
  colnames(stack_bars) <- c("Treshold", "D_type")
  
  barplot <- ggplot(data = stack_bars) + 
    geom_bar(mapping = aes(x = D_type, fill = Treshold), position = "fill") +
    scale_fill_manual("Legend", values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"))+
    theme_classic()
  
  assign("path_exp", paste0("Figures/Barplot_V4D", i, ".pdf"))
  
  pdf(path_exp)
  plot(barplot)
  dev.off()
  
}


## Doing the statistics comparing post 1st vacination and post second

t.test(as.logical(a_1), as.logical(b_1))
t.test(as.logical(a_2), as.logical(b_2))
t.test(as.logical(a_3), as.logical(b_3))
t.test(as.logical(a_4), as.logical(b_4))



# Finding this treshold distribution but separated by serostatus baseline


for(i in 1:length(tresholds_V4)){
  
  SB_Treshold <- data.frame(as.character(tsne_demographics$`Serostatus Baseline`), tresholds_V4[[i]])
  colnames(SB_Treshold) <- c("baseline", "Treshold")
  
  assign("path_exp", paste0("Figures/Threshold_BLNeg_V4D", i, ".pdf"))
  
  pdf(path_exp)
  pie(table(SB_Treshold[SB_Treshold$baseline == 0, ]), col = brewer.set1(4)[3:4], labels = c("Under Treshold", "Above Treshold"))
  dev.off()
  
  print(table(SB_Treshold[SB_Treshold$baseline == 0, ]))
  
  assign("path_exp", paste0("Figures/Threshold_BLPos_V4D", i, ".pdf"))
  
  pdf(path_exp)
  pie(table(SB_Treshold[SB_Treshold$baseline == 1, ]), col = brewer.set1(4)[3:4], labels = c("Under Treshold", "Above Treshold"))
  dev.off()
  
  print(table(SB_Treshold[SB_Treshold$baseline == 1, ]))
}



# Out of the non-responders for the serotype 1, how many where seropositives?



for(i in 1:length(tresholds_V4)){
  
  SB_Treshold <- data.frame(as.character(tsne_demographics$`Serostatus Baseline`), tresholds_V4[[i]])
  colnames(SB_Treshold) <- c("baseline", "Treshold")
  
  
  barplot <- ggplot(data = SB_Treshold) + 
    geom_bar(mapping = aes(x = Treshold, fill = baseline), position = "fill") +
    scale_fill_manual("Legend", values = c("1" = "#4DAF4A", "0" = "#984EA3"))+
    theme_classic()
  
  assign("path_exp", paste0("Figures/Immunization_SB_V4D", i, ".pdf"))
  
  pdf(path_exp)
  plot(barplot)
  dev.off()
  
}





# Doing the statistics of the comparison between the # of patients that reached the 
#threshold of DENV naive and DENV exposed

#DENV1

SB_Treshold_D1 <- data.frame(as.character(tsne_demographics$`Serostatus Baseline`), tresholds_V4[[1]])
colnames(SB_Treshold_D1) <- c("baseline", "Treshold")


SB_Treshold_D1_AT <- SB_Treshold_D1[SB_Treshold_D1$Treshold == T, ]
SB_Treshold_D1_UT <- SB_Treshold_D1[SB_Treshold_D1$Treshold == F, ]

SB_Treshold_D1_AT <- SB_Treshold_D1_AT$baseline
SB_Treshold_D1_UT <- SB_Treshold_D1_UT$baseline

for(i in 1:length(SB_Treshold_D1_AT)){
  if(SB_Treshold_D1_AT[i] == 1){
    SB_Treshold_D1_AT[i] <- T
  }else if(SB_Treshold_D1_AT[i] == 0){
    SB_Treshold_D1_AT[i] <- F
  }
}

for(i in 1:length(SB_Treshold_D1_UT)){
  if(SB_Treshold_D1_UT[i] == 1){
    SB_Treshold_D1_UT[i] <- T
  }else if(SB_Treshold_D1_UT[i] == 0){
    SB_Treshold_D1_UT[i] <- F
  }
}



t.test(as.logical(SB_Treshold_D1_AT), as.logical(SB_Treshold_D1_UT))

#DENV2


SB_Treshold_D2 <- data.frame(as.character(tsne_demographics$`Serostatus Baseline`), tresholds_V4[[2]])
colnames(SB_Treshold_D2) <- c("baseline", "Treshold")


SB_Treshold_D2_AT <- SB_Treshold_D2[SB_Treshold_D2$Treshold == T, ]
SB_Treshold_D2_UT <- SB_Treshold_D2[SB_Treshold_D2$Treshold == F, ]

SB_Treshold_D2_AT <- SB_Treshold_D2_AT$baseline
SB_Treshold_D2_UT <- SB_Treshold_D2_UT$baseline

for(i in 1:length(SB_Treshold_D2_AT)){
  if(SB_Treshold_D2_AT[i] == 1){
    SB_Treshold_D2_AT[i] <- T
  }else if(SB_Treshold_D2_AT[i] == 0){
    SB_Treshold_D2_AT[i] <- F
  }
}

for(i in 1:length(SB_Treshold_D2_UT)){
  if(SB_Treshold_D2_UT[i] == 1){
    SB_Treshold_D2_UT[i] <- T
  }else if(SB_Treshold_D2_UT[i] == 0){
    SB_Treshold_D2_UT[i] <- F
  }
}



t.test(as.logical(SB_Treshold_D2_AT), as.logical(SB_Treshold_D2_UT))





#DENV3


SB_Treshold_D3 <- data.frame(as.character(tsne_demographics$`Serostatus Baseline`), tresholds_V4[[3]])
colnames(SB_Treshold_D3) <- c("baseline", "Treshold")


SB_Treshold_D3_AT <- SB_Treshold_D3[SB_Treshold_D3$Treshold == T, ]
SB_Treshold_D3_UT <- SB_Treshold_D3[SB_Treshold_D3$Treshold == F, ]

SB_Treshold_D3_AT <- SB_Treshold_D3_AT$baseline
SB_Treshold_D3_UT <- SB_Treshold_D3_UT$baseline

for(i in 1:length(SB_Treshold_D3_AT)){
  if(SB_Treshold_D3_AT[i] == 1){
    SB_Treshold_D3_AT[i] <- T
  }else if(SB_Treshold_D3_AT[i] == 0){
    SB_Treshold_D3_AT[i] <- F
  }
}

for(i in 1:length(SB_Treshold_D3_UT)){
  if(SB_Treshold_D3_UT[i] == 1){
    SB_Treshold_D3_UT[i] <- T
  }else if(SB_Treshold_D3_UT[i] == 0){
    SB_Treshold_D3_UT[i] <- F
  }
}



t.test(as.logical(SB_Treshold_D3_AT), as.logical(SB_Treshold_D3_UT))




#DENV4


SB_Treshold_D4 <- data.frame(as.character(tsne_demographics$`Serostatus Baseline`), tresholds_V4[[4]])
colnames(SB_Treshold_D4) <- c("baseline", "Treshold")


SB_Treshold_D4_AT <- SB_Treshold_D4[SB_Treshold_D4$Treshold == T, ]
SB_Treshold_D4_UT <- SB_Treshold_D4[SB_Treshold_D4$Treshold == F, ]

SB_Treshold_D4_AT <- SB_Treshold_D4_AT$baseline
SB_Treshold_D4_UT <- SB_Treshold_D4_UT$baseline

for(i in 1:length(SB_Treshold_D4_AT)){
  if(SB_Treshold_D4_AT[i] == 1){
    SB_Treshold_D4_AT[i] <- T
  }else if(SB_Treshold_D4_AT[i] == 0){
    SB_Treshold_D4_AT[i] <- F
  }
}

for(i in 1:length(SB_Treshold_D4_UT)){
  if(SB_Treshold_D4_UT[i] == 1){
    SB_Treshold_D4_UT[i] <- T
  }else if(SB_Treshold_D4_UT[i] == 0){
    SB_Treshold_D4_UT[i] <- F
  }
}



t.test(as.logical(SB_Treshold_D4_AT), as.logical(SB_Treshold_D4_UT))





### Finding How many patients became immunized against 1, 2, 3, 4 or 0 serotypes



responses <- cbind(as.data.frame(b_1), as.data.frame(b_2))
responses <- cbind(responses, as.data.frame(b_3))
responses <- cbind(responses, as.data.frame(b_4))

groups <- list(c(0), c(0), c(0), c(0),c(0))

for(i in 1:nrow(responses)){
  if(sum(responses[i,]) == 1){
    groups[[1]] <- groups[[1]] + 1
  } else if(sum(responses[i,]) == 2){
    groups[[2]] <- groups[[2]]+1
  } else if(sum(responses[i,]) == 3){
    groups[[3]] <- groups[[3]]+1
  } else if(sum(responses[i,]) == 4){
    groups[[4]] <- groups[[4]]+1
  } else if(sum(responses[i,]) == 0){
    groups[[5]] <- groups[[5]]+1
  }
}


immunized <- c(groups[[1]], groups[[2]], groups[[3]], groups[[4]], groups[[5]])



pdf("Figures/Pie_Immunized_patients.pdf")
pie(immunized, col = brewer.accent(5))
dev.off()

#Creating a column that says if the patient responded and to how many do they respond






re <- c()
for(i in 1:nrow(responses)){
  if(sum(responses[i,]) == 0){
    re[i] <- 0
  }else if(sum(responses[i,]) == 1){
    re[i] <- 1
  }else if(sum(responses[i,]) == 2){
    re[i] <- 2
  }else if(sum(responses[i,]) == 3){
    re[i] <- 3
  }else if(sum(responses[i,]) == 4){
    re[i] <- 4
  }
}

names <- tsne_demographics$names
Threshold <- data.frame(names, re)


tsne_demographics_3 <- merge(tsne_demographics, Threshold, by = "names")


#saveRDS(tsne_demographics_3,"RDS_objects/tsne_demographics_3.rds")
#tsne_demographics_3 <- readRDS("RDS_objects/tsne_demographics_3.rds")



responses$SB <- tsne_demographics_3$`Serostatus Baseline`
responses$X <- rep("a", nrow(responses))

responses_soropos <- responses[responses$SB == 1,] 

responses_neg <- responses[responses$SB == 0,] 

sum(responses_soropos$b_1)/nrow(responses_soropos)
sum(responses_neg$b_1)/nrow(responses_neg)


sum(responses_soropos$b_2)/nrow(responses_soropos)
sum(responses_neg$b_2)/nrow(responses_neg)

sum(responses_soropos$b_3)/nrow(responses_soropos)
sum(responses_neg$b_3)/nrow(responses_neg)

sum(responses_soropos$b_4)/nrow(responses_soropos)
sum(responses_neg$b_4)/nrow(responses_neg)


copy_responses <- responses

for(i in 1:nrow(responses)){
  for(j in 1:ncol(responses)){
    if(responses[i, j] == T){
      responses[i, j] <- "AT"
    } else{
      responses[i, j] <- "UT"
    }
  }
}

responses$SB <- tsne_demographics_3$`Serostatus Baseline`
responses$X <- rep("a", nrow(responses))

responses_soropos <- responses[responses$SB == 1,] 

responses_neg <- responses[responses$SB == 0,] 






pdf("Figures/imprinting_1.pdf")
ggplot(data = responses_soropos) + 
  geom_bar(mapping = aes(x = X, fill = b_1), position = "fill") +
  scale_fill_manual("Legend", values = c("AT" = "#A74E56", "UT" = "#BEDFDA"))+
  theme_classic()
dev.off()

pdf("Figures/imprinting_2.pdf")
ggplot(data = responses_soropos) + 
  geom_bar(mapping = aes(x = X, fill = b_2), position = "fill") +
  scale_fill_manual("Legend", values = c("AT" = "#A74E56", "UT" = "#BEDFDA"))+
  theme_classic()
dev.off()

pdf("Figures/imprinting_3.pdf")
ggplot(data = responses_soropos) + 
  geom_bar(mapping = aes(x = X, fill = b_3), position = "fill") +
  scale_fill_manual("Legend", values = c("AT" = "#A74E56", "UT" = "#BEDFDA"))+
  theme_classic()
dev.off()

pdf("Figures/imprinting_4.pdf")
ggplot(data = responses_soropos) + 
  geom_bar(mapping = aes(x = X, fill = b_4), position = "fill") +
  scale_fill_manual("Legend", values = c("AT" = "#A74E56", "UT" = "#BEDFDA"))+
  theme_classic()
dev.off()

pdf("Figures/imprinting_5.pdf")
ggplot(data = responses_neg) + 
  geom_bar(mapping = aes(x = X, fill = b_1), position = "fill") +
  scale_fill_manual("Legend", values = c("AT" = "#A74E56", "UT" = "#BEDFDA"))+
  theme_classic()
dev.off()

pdf("Figures/imprinting_6.pdf")
ggplot(data = responses_neg) + 
  geom_bar(mapping = aes(x = X, fill = b_2), position = "fill") +
  scale_fill_manual("Legend", values = c("AT" = "#A74E56", "UT" = "#BEDFDA"))+
  theme_classic()
dev.off()

pdf("Figures/imprinting_7.pdf")
ggplot(data = responses_neg) + 
  geom_bar(mapping = aes(x = X, fill = b_3), position = "fill") +
  scale_fill_manual("Legend", values = c("AT" = "#A74E56", "UT" = "#BEDFDA"))+
  theme_classic()
dev.off()

pdf("Figures/imprinting_8.pdf")
ggplot(data = responses_neg) + 
  geom_bar(mapping = aes(x = X, fill = b_4), position = "fill") +
  scale_fill_manual("Legend", values = c("AT" = "#A74E56", "UT" = "#BEDFDA"))+
  theme_classic()
dev.off()



tsne_demographics_3$re <- as.factor(tsne_demographics_3$re)
barplot <- ggplot(data = tsne_demographics_3) + 
  geom_bar(mapping = aes(x = hclust_5, fill = re), position = "fill") +
  theme_classic()


pdf("Figures/TresholdPerCluster.pdf")
barplot
dev.off()



#### Which are the most commom combinations of 3 serotypes immunization



index_threeOrMore <- c()
for(i in 1:nrow(copy_responses)){
  if(sum(copy_responses[i,1:4]) >= 3 ){
    index_threeOrMore <- c(index_threeOrMore, i)
  }
}

copy_responses$baseline <- tsne_demographics$`Serostatus Baseline`



threeOrMore <- copy_responses[index_threeOrMore, ]
table(threeOrMore$baseline)

pdf("Figures/Baseline_tripleImmunization.pdf")
pie(c(17, 30))
dev.off()
17/47
30/47



threeOrMore_pos <- threeOrMore[threeOrMore$baseline == 1, ]
table(threeOrMore_pos$b_4)

threeOrMore_neg <- threeOrMore[threeOrMore$baseline == 0, ]
table(threeOrMore_neg$b_4)



#### Which is the most commomn solo immunization 

copy_responses$baseline <- NULL

only_one <- c()
for(i in 1:nrow(copy_responses)){
  if(sum(copy_responses[i,1:4]) == 1 ){
    only_one <- c(only_one, i)
  }
}

copy_responses$baseline <- tsne_demographics$`Serostatus Baseline`


only_one_df <- copy_responses[only_one, ]

table(only_one_df$baseline)


only_one_type <- c()
for(i in 1:nrow(only_one_df)){
  if(only_one_df$b_1[i] == T ){
    only_one_type <- c(only_one_type, "A")
  } else if(only_one_df$b_2[i] == T ){
    only_one_type <- c(only_one_type, "B")
  } else if(only_one_df$b_3[i] == T ){
    only_one_type <- c(only_one_type, "C")
  } else if(only_one_df$b_4[i] == T ){
    only_one_type <- c(only_one_type, "D")
  }
}


table(only_one_type)
only_one_df$which <-only_one_type


pos <- subset(only_one_df, subset = baseline == 1)
neg <- subset(only_one_df, subset = baseline == 0)

table(pos$which)
which_pos <- c(2, 0, 1)


table(neg$which)
which_neg <- c(1, 11, 2)

df_treemap <- data.frame(Baseline = c(rep("Positive", 3), rep("Negative", 3)), 
                         Serotype = c("DENV1", "DENV2", "DENV3", "DENV1", "DENV2", "DENV3"), 
                         Values = c(which_pos, which_neg))


pdf("Figures/One_Positive.pdf")
treemap(df_treemap, index = c("Baseline", "Serotype"),
        vSize = "Values",
        border.col = c( "white"), # Color of borders of groups, of subgroups, of subsubgroups ....
        palette = "Dark2",
        border.lwds = c(7, 2) # Width of colors
)

dev.off()





a <- data.frame(only_one_type, rep("x", length(only_one_type)))
colnames(a) <- c("type", "x")

barplot <- ggplot(data = a) + 
  geom_bar(mapping = aes(x = x, fill = type), position = "fill") +
  scale_fill_manual("Legend", values = c("A"="#A6CEE3", "B"="#B2DF8A" , "C"="#FB9A99"))+
  theme_classic()











####### ---------------------- Randon Forest to study feature importace





demographic_fr <- na.omit(tsne_demographics[,8:18 ])
demographic_fr$ID <-  NULL
demographic_fr$IDs <-  NULL
demographic_fr$hclust_5 <- as.factor(demographic_fr$hclust_5)
demographic_fr$`Serostatus Baseline BR` <- NULL
demographic_fr$`Age Decaded` <- NULL


head(demographic_fr)
# Changing the names because the randonforest package can't deal with double names 

colnames(demographic_fr) <- c("Cluster", "DD", "VFA", "Age", "SB", "Sex", "Age_2")
demographic_fr$Age_2 <- NULL
head(demographic_fr)

#Splitting the data between train and test sets 

source("Functions/trainTestSplit.R")

model_data <- trainTestSplit(demographic_fr, 84, 10)
train_set <- model_data[[1]]
test_set <- model_data[[2]]

# Creating and optimizing the model

source("Functions/optimize_RandForest.R")

ntree_vector <- c(5, 10, 15, 20, 25, 30, 40, 50)

rf_list <- optimize_RandForest("Cluster", train_set, test_set, ntree_vector)

final_model <- rf_list[["model"]]

varPlot <- rf_list[["FI_plot"]]


pdf("Figures/RF_85Acu.pdf")
varPlot
dev.off()






