###################################
#### Author: Raúl Mejía
#### The aim of this program is to explore your expression matrix, taking a glimpse of the normalization & batch effects  
###################################
# one column of your annotation file should have the name "group" and other "Unique_ID"
# the column names of your expression matrix should match with the rows entries of the "group" column in your annotation df (in that exact order)
# This script calls the function PCA_box_density_plots_group_Treatment_Cell_line that plot the group Treatment and Cell_line columns from your annotation file (those columns should be in your annotation file)
###################################
#### 0) loading and/or installing required libraries
###################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("pheatmap")) {
  BiocManager::install("pheatmap", ask =FALSE)
  library("pheatmap")
}
if (!require("ggplot2")) {
  BiocManager::install("ggplot2", ask =FALSE)
  library("ggplot2")
}
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("ggfortify")) {
  install.packages("ggfortify", ask =FALSE)
  library("ggfortify")
}
if (!require("sva")) {
  BiocManager::install("sva", ask =FALSE)
  library("sva")
}
if (!require("limma")) {
  BiocManager::install("limma", ask =FALSE)
  library("limma")
}
if (!require("smooth")) {
  install.packages("smooth", ask =FALSE)
  library("smooth")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", ask =FALSE)
  library("RColorBrewer")
}
if (!require("plotrix")) {
  install.packages("plotrix", ask =FALSE)
  library("plotrix")
}
if (!require("datasets")) {
  install.packages("datasets", ask =FALSE)
  library("datasets")
}
if (!require("reshape2")) {
  install.packages("reshape2", ask =FALSE)
  library("reshape2")
}
if (!require("ggridges")) {
  install.packages("ggridges", ask =FALSE)
  library("ggridges")
}
if (!require("cowplot")) {
  install.packages("cowplot", ask =FALSE)
  library("cowplot")
}
if (!require("preprocessCore")) {
  BiocManager::install("preprocessCore", ask =FALSE)
  library("preprocessCore")
}
if (!require("affy")) {
  BiocManager::install("affy", ask =FALSE)
  library("affy")
}
if (!require("oligo")) {
  BiocManager::install("oligo", ask =FALSE)
  library("oligo")
}
if (!require("Rtsne")) {
  BiocManager::install("Rtsne", ask =FALSE)
  library("Rtsne")
}
#if (!require("M3C")) {
#  BiocManager::install("M3C", ask =FALSE)
#  library("M3C")
#}
#if (!require("tidyverse")) {
#  BiocManager::install("tidyverse", ask =FALSE)
#  library("tidyverse")
#}

###################################
#### Data given by the user
###################################
myargs <- commandArgs(trailingOnly = TRUE)
path_to_your_matrix <- myargs[1]
# path_to_your_matrix <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29_to_work_colsymb-n-gene-deleted.tsv"
# path_to_your_matrix <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29-switched-to_work_colsymb-n-gene-deleted_original_order.tsv"

path_to_your_annotation_file <- myargs[2]
# path_to_your_annotation_file <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/annotation_lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_09_Rformat.tsv"

Code_path <- myargs[3] # Path where are the rest of your scripts
# Code_path <- "/media/rmejia/mountme88/code/Processing_an_RNAexpression_matrix_from_hyb_n_scanning"  

path_Results_directory <- myargs[4]
# path_Results_directory <-"/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Results_from_the_exploratory_analysis"

data_label <- myargs[5]
# data_label<- "lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29_to_work_colsymb-n-gene-deleted"
# data_label<- "switched_org_orderlipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29_to_work_colsymb-n-gene-deleted"

colname_4_intra_batch_normalization <- myargs[6] # please don't use spaces or parenthesis/brackets in the names of your columns
# colname_4_intra_batch_normalization <- "group" # the name of your column to correct


###################################
#### Creating your result folders if they doesn't exist yet
###################################
dir.create(path_Results_directory , recursive = TRUE)

###################################
#### Normalize your paths
###################################
Code_path <- normalizePath(Code_path)
path_Results_directory  <- normalizePath( path_Results_directory  )

###################################
#### Reading the annotation table and the table that cointains the expression data
###################################
Raw_expmat <- read.table( path_to_your_matrix  , sep = "\t", header = TRUE)
annot <- read.table( path_to_your_annotation_file , sep = "\t", header = TRUE )

#####################
# Annotation object for plotting pcas
####################
annot_4_plotting_pca <- annot
annot_4_plotting_pca[ , "group" ] <- as.factor( annot_4_plotting_pca[ , "group" ] )
colnames(Raw_expmat) ==  annot_4_plotting_pca$Unique_ID

# loading the function to melt (reshape) the data to preparation for ggplot2 functions
source( paste0( Code_path,"/libraries/","matrix_N_annotdf_2_melteddf.R") )
meltedrawdata <- matrix_N_annotdf_2_melteddf( Raw_expmat , annot )

########################################
########
########    0. Exploratory
########
########################################
#########
### Visualize your the Raw data
########
source(paste0( Code_path,"/libraries/","PCA_box_density_plots.R") )
PCA_box_density_plots_group_Treatment_Cell_line(  paste0( path_Results_directory,"/Exploratory" )  ,
                        Raw_expmat ,  annot_4_plotting_pca , meltedrawdata , paste0( data_label, "Data_as_given" ))

########################################
########
########    1. Preprocessing
########
########################################
##################
## log 2 transformation
##################
expmat_log2 <- log2( Raw_expmat +1)
melted_expmat_log2 <- matrix_N_annotdf_2_melteddf(expmat_log2 , annot )

# Visualize the data log2 transformed data
PCA_box_density_plots_group_Treatment_Cell_line(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        expmat_log2 ,  annot_4_plotting_pca , melted_expmat_log2 , paste0( data_label, "_Log2" ))

write.table(expmat_log2, file= paste0(path_to_your_matrix ,"_log2.tsv") , sep="\t", row.names = TRUE, col.names = TRUE )