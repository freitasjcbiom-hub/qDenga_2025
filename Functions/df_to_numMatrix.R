################################
# Function to transform a df with num as character to numeric,

df_to_numMatrix <- function(df) {
  cols <- colnames(df)
  rows <- rownames(df)
  matrix_data <- matrix(0, nrow = nrow(df), ncol = ncol(df))
  for(i in 1:ncol(df)){
    matrix_data[, i] <- c(as.numeric(unlist(df[,i])))
  }
  colnames(matrix_data) <- cols
  rownames(matrix_data) <- rows
  return(matrix_data)
}