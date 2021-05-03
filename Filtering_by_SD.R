# This program eliminates the features that doesn´t fulfill the criteria
# The criteria is to have a SD / SD mx < chosen parameter
# The SD is calculated by condition
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

SD_ratio <- myargs[5]
# SD_ratio <- 4  

data_label <- myargs[6]
# data_label<- SD_stable_features_per_condition

#############
## Reading the data
#############
dir.create(path_Results_directory, recursive = TRUE)
path_Results_directory <- normalizePath(path_Results_directory)
mymatrix <- read.table( path_to_your_matrix , sep = "\t", header = TRUE )
annotdf <- read.table(path_to_your_annotation_file , sep = "\t", header = TRUE )

if(all(colnames(mymatrix) ==  annotdf$Unique_ID) == TRUE ){ 
  print("your annotation and colnames match")
  }
if(all(colnames(mymatrix) ==  annotdf$Unique_ID) != TRUE ){ 
  print("ERROR: your annotation and colnames DON´T match")
  break()
}





split( )