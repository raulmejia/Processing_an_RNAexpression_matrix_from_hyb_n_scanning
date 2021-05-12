# Input a numeric matrix, an annotation file, and a row of the matrix that we want to generate the box plot
# One column in your annot file should be called Unique_ID
# The name "values" is reserved for the columns 
# This script has the following restriction: expmatrix´s colnames and Unique_ID from your annotation file match and in the same order,
#     Maybe you can edit it to ask only for the same names but not in the same order
#
# Example of use:
# Rscript Box_plot_element_from_this_matrix_under_this_annotations.R /path/your_matrix.tsv /path/your_annotation_file.tsv /path/to/code /folder/path/for/your/results Exact_Name_of_your_row Title_for_your_plot Name_of_the_source_matrix_to_stamp_it_in_the_pdf
#
# Note:
# I usually receive this message : `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
# That means that a default value for the number of bins was choosen for your data and that maybe is not the best for your data.
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
if (!require("tidyverse")) {
  BiocManager::install("tidyverse", ask =FALSE)
  library("tidyverse")
}
if (!require("hrbrthemes")) {
  BiocManager::install("hrbrthemes", ask =FALSE)
  library("hrbrthemes")
}
if (!require("viridis")) {
  BiocManager::install("viridis", ask =FALSE)
  library("viridis")
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
# row_to_plot <- "ENSG00000132518"

data_label <- myargs[6]
# data_label<- "GUCY2D"
# data_label<- "switched_org_orderlipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29_to_work_colsymb-n-gene-deleted"

label_matrix_of_origin <- myargs[7]
# label_matrix_of_origin <- "lipidosis_RNA_16_STAR_fC_edgeR_matrix_switched_log2"

#############
## Reading the data
#############
dir.create(path_Results_directory, recursive = TRUE)
path_Results_directory <- normalizePath(path_Results_directory)
mymatrix <- read.table( path_to_your_matrix , sep = "\t", header = TRUE )
annotdf <- read.table(path_to_your_annotation_file , sep = "\t", header = TRUE )

values_toplot <- mymatrix[grep(row_to_plot , rownames(mymatrix) ), ]

# Check if the annotation file is appropiate
if( all(annotdf[,"Unique_ID"] == names(values_toplot)) == TRUE){ 
  print("expmatrix´s colnames and Unique_ID from your annotation file match and in the same order")
}

if( all(annotdf[,"Unique_ID"] != names(values_toplot)) == TRUE){ 
  print("ERROR! the names doesn't match")
}

rownames(annotdf) <- names( values_toplot )
object_to_plot <- data.frame( as.numeric(values_toplot) , annotdf )
colnames(object_to_plot)[1] <- "values"

pdf_file_path <-paste0( path_Results_directory , "/" ,"Plotted-row_",row_to_plot,"_Labeled-as_", data_label ,"_From-the-matrix__",label_matrix_of_origin , ".pdf" )

pdf(file= pdf_file_path )
ggplot(object_to_plot, aes(x=group, y=values, fill=Treatment)) + 
  scale_fill_viridis(discrete = TRUE, alpha=0.4)+
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(data_label) + theme_grey(base_size = 22)
dev.off()






