###################################################
## Function to optimize the randonforest parameter#


optimize_RandForest <- function(category, trainData, testData, n_tree_vector){
  set.seed(123) # setting a seed for reproducibility
  library(randomForest) # loading the required package
  
  formula <- as.formula(paste(category, "~ .")) # Dynamically create the formula using the target column name

  accuracy_vector <- c() # A vector to store the accuracies
  
  
  for(i in 1:length(n_tree_vector)){ # Creating the model and calculating the accuracy for each one of the possible ntrees
    
    model <- randomForest(formula, data = trainData, ntree = n_tree_vector[i] ) #Creating the model
    
    test_labels <- testData[[category]]
    test_length <- length(test_labels)
    
    predictions <- predict(model, testData) #Making prediction
    
    accuracy_val <- sum(predictions == test_labels) / test_length #Checking the accuracy
    
    accuracy_vector[i] <- accuracy_val
    
  }
  
  # Selecting the best accuracy
  maximum_acuracy <- max(accuracy_vector)
  idx_max_acc <- which(accuracy_vector == max(accuracy_vector))
  
  # Printing which was the best accuracy
  print(paste0("The best accuracy is ", maximum_acuracy*100, "%, generated when using ", n_tree_vector[idx_max_acc], " trees."))
  
  
  # Generating the model with the optimum tree number
  model <- randomForest(formula, data = trainData, ntree = n_tree_vector[i] )
  features_importance <- importance(model)
  features_importance_plot <- varImpPlot(model)
  
  # Creating a list to return the results
  results <- list("model" = model, "feat_import" = features_importance, "FI_plot" = features_importance_plot)
  
  # Returning
  return(results)
}








