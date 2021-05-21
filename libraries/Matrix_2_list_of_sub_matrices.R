bye1stcol <- function(df){
  print(dim(df))
  df<- df[,-1]
  print(dim(df))
  return(df)
} 
putcolnmaes <- function(df){
  colnames(df) <- colnames(mymatrix)
  return(df)
}
makeme_num_mat <-function(df){
  df <- as.matrix(df)
  mode(df) <- "numeric"
  return(df)
}


Matrix_2_list_of_sub_matrices <- function( factor_2_split_amatrix, expmat){
  amatrix <- t(expmat)
  matrix4intrabatchnorm <- cbind(as.character( factor_2_split_amatrix ),as.data.frame( amatrix ))
  colnames(matrix4intrabatchnorm)[1] <- "batch"
  matrix4intrabatchnorm[,"batch"] <- as.factor( matrix4intrabatchnorm[,"batch"] )
    list_splited_table <- split( matrix4intrabatchnorm, matrix4intrabatchnorm$batch )
    list_splited_table_no1stcol <- lapply(list_splited_table, bye1stcol)
    list_splited_table_t <- lapply(list_splited_table_no1stcol, t )
  
  
  list_splited_nummat <- lapply(list_splited_table_t, makeme_num_mat)
  
  return(list_splited_nummat)
}
