###########################################
###########################################
## Analysis of eficaccy - Takeda Vaccine ##
###########################################
################ AIM DATA #################
################ T Cells ##################
###########################################



library(readxl)
library(PCAtools)
library(Rtsne)
library(factoextra)
library(patchwork)
library(alluvial)
library(dplyr)
library(ggalluvial)





#---------------------- Correlation between aim and neutralization

# Loading new AIM data 

data_aim <- read_excel("Data/DENV_T_Cell_AIM_Wide.xlsx", sheet = 4)



# Creating a vector to say if that patient is above or bellow the threshold


cutoff_CD4 <- c()
for(i in 1:nrow(data_aim)){
  if(data_aim$CD4_AIM_Dif[i] >= 0.138){
    cutoff_CD4[i] <- 1 # Recieves 1 if it is above the cutoff
  } else{
    cutoff_CD4[i] <- 0 # Recieves 0 if it is above the cutoff
  }
}




cutoff_CD8 <- c()
for(i in 1:nrow(data_aim)){
  if(data_aim$CD8_AIM_Dif[i] >= 0.362){
    cutoff_CD8[i] <- 1 # Recieves 1 if it is above the cutoff
  } else{
    cutoff_CD8[i] <- 0 # Recieves 0 if it is above the cutoff
  }
}





# Creating a new vector that combine the information of the two vectors above, to create an alluvial plot



type_of_response <- c()

for(i in 1:length(cutoff_CD4)){
  if(cutoff_CD4[i] + cutoff_CD8[i] == 2){
    type_of_response[i] <- "TCD4 and TCD8"
  } else if(cutoff_CD4[i] == 1){
    type_of_response[i] <- "TCD4"
  } else if(cutoff_CD8[i] == 1){
    type_of_response[i] <- "TCD8"
  } else{
    type_of_response[i] <- "No Response"
  }
}



# Combining it to the patient ID in a new data frame


ID <- data_aim$sample_ID
type_of_response <- as.factor(type_of_response)


aim_cutoff <- data.frame(ID, type_of_response)




# Loading the demographics data for integration

#tsne_demographics_3 <- readRDS("RDS_objects/tsne_demographics_3.rds")

# Merging it

tsne_demographics_4 <- merge(tsne_demographics_3, aim_cutoff, by = "ID")

# Creating the alluvial


library(dplyr)

data_freq <- tsne_demographics_4 %>%
  group_by(re, type_of_response) %>%
  summarise(Freq = n(), .groups = "drop")


alluvial(
  data_freq[, 1:2],     # the categorical variables
  freq = data_freq$Freq, # flow weight
  border = NA,
  alpha = 0.8)

percent_response <- c()
for(i in c(1, 5, 9, 13, 17)){
  total <-  sum(data_freq$Freq[i],
                data_freq$Freq[i+1],
                data_freq$Freq[i+2],
                data_freq$Freq[i+3])
  percent_response <- c(percent_response, data_freq$Freq[i]/total)
  percent_response <- c(percent_response, data_freq$Freq[i+1]/total)
  percent_response <- c(percent_response, data_freq$Freq[i+2]/total)
  percent_response <- c(percent_response, data_freq$Freq[i+3]/total)
}


alluvial(
  data_freq[, 1:2],     
  freq = percent_response, 
  border = NA,
  col = brewer.set1(4),
  alpha = 0.8)


pdf("Figures/Alluvial_tResponse.pdf", width = 12)
alluvial(
  data_freq[, 1:2],     # the categorical variables
  freq = percent_response, # flow weight
  border = NA,
  col = brewer.set1(4),
  alpha = 0.8)
dev.off()

data_freq$percentage <- percent_response

write.csv(data_freq, "Data/AIM_and_Neutralization.csv")




data_freq <- tsne_demographics_4 %>%
  group_by(`Serostatus Baseline`, type_of_response) %>%
  summarise(Freq = n(), .groups = "drop")


alluvial(
  data_freq[, 1:2],     # the categorical variables
  freq = data_freq$Freq, # flow weight
  border = NA,
  col = brewer.set1(10)[3:6],
  alpha = 0.8)




pdf("Figures/Alluvial_tResponse_baseline.pdf", width = 12)
alluvial(
  data_freq[, 1:2],     # the categorical variables
  freq = data_freq$Freq, # flow weight
  border = NA,
  col = brewer.set1(10)[3:6],
  alpha = 0.8)
dev.off()


# Changing the graph as requested:

# Doing Tcell responders x Non responders


aim_cutoff <- data.frame(ID, type_of_response)

binary_t_response <- c()
for(i in 1:length(type_of_response)){
  if(type_of_response[i] == "No Response"){
    binary_t_response[i] <- 0
  } else{
    binary_t_response[i] <- 1
  }
}


binary_t_response <- as.factor(binary_t_response)

aim_cutoff_2 <- data.frame(ID, binary_t_response)

tsne_demographics_4 <- merge(tsne_demographics_3, aim_cutoff_2, by = "ID")



data_freq <- tsne_demographics_4 %>%
  group_by(re, binary_t_response) %>%
  summarise(Freq = n(), .groups = "drop")



pdf("Figures/Alluvial_tResponse_binary.pdf")
alluvial(
  data_freq[, 1:2],     # the categorical variables
  freq = data_freq$Freq, # flow weight
  col = brewer.set1(10)[3:4],
  border = NA,
  alpha = 0.8)
dev.off()


percent_response_bin <- c()
for(i in seq(1, nrow(data_freq), 2)){
  total <-  sum(data_freq$Freq[i],
                data_freq$Freq[i+1])
  percent_response_bin <- c(percent_response_bin, data_freq$Freq[i]/total)
  percent_response_bin <- c(percent_response_bin, data_freq$Freq[i+1]/total)
}



pdf("Figures/Alluvial_tResponse_binary_normalized.pdf")
alluvial(
  data_freq[, 1:2],     # the categorical variables
  freq = percent_response_bin, # flow weight
  col = brewer.set1(10)[3:4],
  border = NA,
  alpha = 0.8)
dev.off()

data_freq$percentage <- percent_response_bin



