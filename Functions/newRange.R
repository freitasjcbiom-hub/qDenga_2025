# Function to delimitate upper and lower boundaries in a df 


newRange <- function(matrix, upper_boundarie, lower_boundarie){
  df <- as.data.frame(matrix)
  
  for(i in 1:nrow(df)){
    for(j in 1:ncol(df)){
      if(df[i, j] > upper_boundarie){
        df[i, j] <-  upper_boundarie
      } else if(df[i, j] < lower_boundarie) {
        df[i, j] <-  lower_boundarie
      }
    }
  }
  new_matrix <- as.matrix(df)
  
  return(new_matrix)
}