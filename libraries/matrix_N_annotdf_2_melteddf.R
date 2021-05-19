# This program receives a matrix (genes as rows and samples as columns) as well as an annotation dataframe
# The annotation data frame should have a column called "Unique ID", ids will be taken from that column.
# The program will retrieve a "Melted data frame" with the "Unique_ID" column as the ids

#onematrix <- mymatrix
#oneannotdf <- annotdf

matrix_N_annotdf_2_melteddf <- function( onematrix, oneannotdf){
  
  onematrix_t <- t(onematrix)
  table_with_uniqID <- cbind( oneannotdf , as.data.frame(onematrix_t) )
  
  #  table_with_uniqID <- cbind( as.character(oneannotdf$Unique_ID) , as.data.frame(onematrix_t) )
  #  colnames(table_with_uniqID)[1] <- c("Unique_ID")
  
  table_melted_MtxAndUniqueID <- melt( data=table_with_uniqID, 
                                       id.vars=intersect( colnames( table_with_uniqID ) , colnames(oneannotdf) ), 
                                       measure.vars = setdiff( colnames( table_with_uniqID ) , colnames(oneannotdf) )   )
  
  #table_melted_MtxAndUniqueID <- melt( data=table_with_uniqID, id.vars="Unique_ID", measure.vars = colnames(table_with_uniqID)[-1]   )
  
  #colnames(table_melted_MtxAndUniqueID)[2] <- "gene"
  
 # group <-                  as.character( rep( oneannotdf$group, length(unique(table_melted_MtxAndUniqueID[,"gene"]))) )
  
  #result_table_melted_MtxAndScanID <- cbind(table_melted_MtxAndUniqueID , group)
  
  return(table_melted_MtxAndUniqueID)
}
