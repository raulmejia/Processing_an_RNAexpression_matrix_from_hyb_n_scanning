# this program splits a matrix according with an annot file
# the annotation file should have a column named "group"
# and the names of the matrix columns should match with the column "Unique_ID" in the annotation data frame 
split_matrix_by_group_column_in_annotdf <- function(matrix_to_split ,  somedf_with_group_column ){
  list_of_sub_dfs<-list()
  for( k in 1:length(unique(somedf_with_group_column$group))){
    positions <- which(somedf_with_group_column$group %in% unique(somedf_with_group_column$group)[k])
    positions 
    colnames_to_extract <- somedf_with_group_column[positions , "Unique_ID"]
    list_of_sub_dfs[[k]] <- matrix_to_split[, colnames_to_extract]  
    
  }
  names( list_of_sub_dfs) <- unique(somedf_with_group_column$group)
  return(list_of_sub_dfs)
}