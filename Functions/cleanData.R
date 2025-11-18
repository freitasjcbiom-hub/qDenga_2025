#Function to clean a dataframe based on the number of NAs at both rows and collumns


cleanData <- function(df, num_MAX_NA_col, num_MAX_NA_row){ 
  
  col_keep <- rep(FALSE, ncol(df))
  for(i in 1:ncol(df)){
    col_keep[i] <- sum(is.na(df[,i])) <= num_MAX_NA_col
  }
  
  df_filt_1 <- df[ ,col_keep]
  
  row_keep <- rep(FALSE, nrow(df_filt_1))
  for(i in 1:nrow(df_filt_1)){
    row_keep[i] <- sum(is.na(df_filt_1[i,])) <= num_MAX_NA_row 
  }
  
  df_filt_2 <- df_filt_1[row_keep, ]
  
  
  results <- list(df_filt_2, table(col_keep),  table(row_keep)) # Return also the table of each filter
  
  return(results)
}
