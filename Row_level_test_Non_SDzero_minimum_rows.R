# This program eliminates the rows that have SDzero in a maximum number of groups 
# maximun number of groups allowed with SDzero (for each row)
# Your colnames and the ID_Unique column from your annotation file should be the same and in the same order
# Your SD will be calculated according the "group" column in your annot file
# This program calls to the function /libraries/","split_matrix_by_group_column_in_annotdf.R
#
# to do:
# Set the default in an automatic way: Rowlevel_Maximum_number_of_groups_with_SDzero <- length(unique(annotdf$group)) - 2
###################################
###################################
#### 0) loading and/or installing required libraries
###################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("ggplot2")) {
  BiocManager::install("ggplot2", ask =FALSE)
  library("ggplot2")
}
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("limma")) {
  install.packages("limma", ask =FALSE)
  library("limma")
}

###################################
#### Data given by the user
###################################
myargs <- commandArgs(trailingOnly = TRUE)
path_to_your_matrix <- myargs[1]
# path_to_your_matrix <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29-switched-to_work_colsymb-n-gene-deleted_original_order.tsv_log2.tsv"

path_to_your_annotation_file <- myargs[2]
# path_to_your_annotation_file <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/annotation_lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_09_Rformat.tsv"

Code_path <- myargs[3] # Path where are the rest of your scripts
# Code_path <- "/media/rmejia/mountme88/code/Processing_an_RNAexpression_matrix_from_hyb_n_scanning"  

path_Results_directory <- myargs[4]
# path_Results_directory <-"/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Box_Plot_particular_Genes"

Rowlevel_Maximum_number_of_groups_with_SDzero <- myargs[5]
# The Default is to the number of groups minus 2 (In other words at least 2 non-SDzero at rowLevel)
# Rowlevel_Maximum_number_of_groups_with_SDzero <- length(unique(annotdf$group)) - 2 # Remember to set the default in an automatic way
# Rowlevel_Maximum_number_of_groups_with_SDzero <- 2

data_label <- myargs[6]
# data_label<- "non_zeroSD_row_level"

#############
## Reading the data
#############
dir.create(path_Results_directory, recursive = TRUE)
path_Results_directory <- normalizePath(path_Results_directory)
Code_path <- normalizePath(Code_path)
mymatrix <- read.table( path_to_your_matrix , sep = "\t", header = TRUE )
annotdf <- read.table(path_to_your_annotation_file , sep = "\t", header = TRUE )

if(all(colnames(mymatrix) ==  annotdf$Unique_ID) == TRUE ){ 
  print("your annotation and colnames match")
}
if(all(colnames(mymatrix) ==  annotdf$Unique_ID) != TRUE ){ 
  print("ERROR: your annotation and colnames DON´T match")
  break()
}

#############
## Working with the data
#############
wholeSDs <- apply(mymatrix,1,sd) # Getting the total SD per rows

source( paste0(Code_path,"/libraries/","split_matrix_by_group_column_in_annotdf.R"))
splited_matrix_bygroup <- split_matrix_by_group_column_in_annotdf( mymatrix, annotdf) # split in submatrices according to group

sd_per_row_per_condition <- lapply( splited_matrix_bygroup , function(x){ apply( x, 1, sd)}) # getting the SD per rows

sd_per_row_per_condition_matrix <- do.call( cbind, sd_per_row_per_condition )
sd_per_row_per_condition_matrix_Nonzero <- sd_per_row_per_condition_matrix
sd_per_row_per_condition_matrix_Nonzero[sd_per_row_per_condition_matrix_Nonzero != 0] <- "Non-zero"
sd_per_row_per_condition_matrix_Nonzero[sd_per_row_per_condition_matrix_Nonzero == 0] <- 1
sd_per_row_per_condition_matrix_Nonzero[sd_per_row_per_condition_matrix_Nonzero == "Non-zero"] <- 0 
mode( sd_per_row_per_condition_matrix_Nonzero ) <- "numeric"

rowSums_sd_per_row_per_condition_matrix_Nonzero <- rowSums( sd_per_row_per_condition_matrix_Nonzero )

to_much_sd_zeros_perrow <- rowSums_sd_per_row_per_condition_matrix_Nonzero > Rowlevel_Maximum_number_of_groups_with_SDzero
positions_with_acceptable_sd_zeros_per_condition <- !to_much_sd_zeros_perrow
mymatrix_SDzero_per_condition_filtered <- mymatrix[positions_with_acceptable_sd_zeros_per_condition,] # successful rows
Matrix_with_The_lost_rows <- mymatrix[to_much_sd_zeros_perrow,] # failed rows

print("number_of_original_rows"); print(dim(mymatrix) )
print("number_of_lost_rows"); print(dim(Matrix_with_The_lost_rows))
print("percentage_of_retained_rows"); print(paste(dim(mymatrix_SDzero_per_condition_filtered  )[1]/dim(mymatrix)[1] *100 ,"%") )

# rows_that_didn´t pass the filter of zeroSD in Minimum group per condition 

failed_path <- paste0( path_to_your_matrix, "_rows_that_failed_the_non_SDzero-test_in_a_minimum_of_",Rowlevel_Maximum_number_of_groups_with_SDzero,"_groups.tsv")
write.table(Matrix_with_The_lost_rows , file=failed_path , sep="\t", col.names=TRUE, row.names= TRUE)

successful_path <- paste0( path_to_your_matrix, "_successful_rows_non_SDzero-test_in_a_minimum_of_",Rowlevel_Maximum_number_of_groups_with_SDzero,"_groups.tsv")
write.table( mymatrix_SDzero_per_condition_filtered , file= successful_path , sep="\t", col.names=TRUE, row.names= TRUE ) 


