# Function to split a given dataset into training and test



trainTestSplit <- function(df, num_train, num_test){
  all_samples <- 1:nrow(df)
  
  train_index <- sample(all_samples, size = num_train, replace = F)
  test_index <- all_samples[!(all_samples %in% train_index)]
  
  train_set <- df[train_index, ]
  test_set <- df[test_index, ]
  
  sets <- list("train_set" <- train_set, "test_set" <- test_set)
  
  return(sets)
}





