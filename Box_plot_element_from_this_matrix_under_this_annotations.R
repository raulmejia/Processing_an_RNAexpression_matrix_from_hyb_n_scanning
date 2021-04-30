# Input a numeric matrix, an annotation file, and a row of the matrix that we want to generate the box plot
# One column in your annot file should be called Unique_ID
# The name "values" is reserved for the columns 
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

row_to_plot <- myargs[5]
# row_to_plot <- "ENSG00000198691"

data_label <- myargs[6]
# data_label<- "ABCA4"
# data_label<- "switched_org_orderlipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29_to_work_colsymb-n-gene-deleted"

#############
## Reading the data
#############
dir.create(path_Results_directory, recursive = TRUE)
path_Results_directory <- normalizePath(path_Results_directory)
mymatrix <- read.table( path_to_your_matrix , sep = "\t", header = TRUE )
annotdf <- read.table(path_to_your_annotation_file , sep = "\t", header = TRUE )

values_toplot <- mymatrix[grep(row_to_plot , rownames(mymatrix) ), ]

if( all(annotdf[,"Unique_ID"] == names(values_toplot)) == TRUE){ 
  print("the names match")
}

if( all(annotdf[,"Unique_ID"] != names(values_toplot)) == TRUE){ 
  print("ERROR! the names doesn't match")
}
rownames(annotdf) <- names( values_toplot )
object_to_plot <- data.frame( as.numeric(values_toplot) , annotdf )
colnames(object_to_plot)[1] <- "values"

ggplot(object_to_plot, aes(x=group, y=values, fill=Treatment)) + 
  geom_boxplot() + ggtitle(data_label) + theme_grey(base_size = 22)



