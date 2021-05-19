# This program eliminates the features that doesn´t fulfill the criteria
# The criteria is to have a SD / SD mx < chosen parameter
# The SD is calculated by condition
# Your colnames and the ID_Unique column from your annotation file should be the same and in the same order
# Your SD will be calculated according the "group" column in your annot file
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

SD_ratio <- myargs[5]
# SD_ratio <- 4  

data_label <- myargs[6]
# data_label<- "SD_stable_features_per_condition"

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
source( paste0(Code_path,"/libraries/","split_matrix_by_group_column_in_annotdf.R"))
splited_matrix_bygroup <- split_matrix_by_group_column_in_annotdf( mymatrix, annotdf)

head(splited_matrix_bygroup[[1]])







wholeSDs <- apply(mymatrix,1,sd)


str(reshaped)
head(reshaped)
dim(reshaped)

58884*4

?do.call
head(wholeSDs)

ratio_absoluteSDs <- wholeSDs/max(wholeSDs)

plot(ratio_absoluteSDs)

apply(mymatrix_splited,  ,  )
plotDensities(t(mymatrix))
plotDensities(mymatrix, legend = FALSE)

?plotDensities

length(wholeSDs)
dim(mymatrix)
head(mymatrix_splited)

test <- c(1,3,2,1,0,20,0,30,20,21,25,19)

sd(test)
sd(test[1:4])
sd(test[5:8])
sd(test[9:12])
