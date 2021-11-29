# This script do some plots over your matrix
# one group of plots is about the distribution of your data ( box and densityplots)
# the other group of plots are for clustering
#   tsnes, pcas, ...

# The structure of your matrix should be 
#          sample1 sample2 ...
# rowname1  0.42    45
# rowname2  1       0
# rowname3  45      0
# NEG_Prob1 0       12
# POS_E     2       1

# example of the annotation: (Note that columns with only numbers give problems)
#
# "Unique_ID"     "group" "Treatment"     "Cell_line"     "Generation"
# "LipChl1"       "LipChl"        "CQ"    "DIP"   "_1"
# "LipChl2"       "LipChl"        "CQ"    "DIP"   "_2"

# one column of your annotation file should have the name "group" (Check if this is still valid) and other "Unique_ID"

# Example of use:
# Rscript /Path/to/this/nice/script/ExpMat_n_Annot_Log2transform_then_Plotting_PCA_Box_Heatmap_density_tsnes.R \
# -m /Path/to/the/expression/matrix.txt \
# -a /Path/2/your/annotation/file \
# -c /Path/to/find/the/libraries \
# -l some_label_for_the_results \
# -g group \
# -p 4 \
# -s 2 \
# -d 5 \
# -w FALSE \
# -o /my/output/folder/2/create

############################## 
## Required libraries
##############################
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
if (!require("preprocessCore")) {
  BiocManager::install("preprocessCore", ask =FALSE)
  library("preprocessCore")
}
if (!require("affy")) {
  BiocManager::install("affy", ask =FALSE)
  library("affy")
}
if (!require("Rtsne")) {
  BiocManager::install("Rtsne", ask =FALSE)
  library("Rtsne")
}

############################## 
## Data given by the user
##############################
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-m", "--matrix", type="character", 
                    help="path to your expression matrix")
parser$add_argument("-a", "--annotation", type="character", 
                    help="path to your annotation file")
parser$add_argument("-c", "--code", type="character", 
                    help="path to your code")
parser$add_argument("-l", "--label", type="character", 
                    help="label to your results")
parser$add_argument("-g", "--maingroups", type="character", 
                    help="the name of your column to correct / make intrabatch normalization")
parser$add_argument("-p", "--perplexity", type="character", 
                    help="Perplexity number for Tsnes")
parser$add_argument("-s", "--pcapointsize", type="character", 
                    help="PCA point size")
parser$add_argument("-w", "--whatapheatmat", type="character", 
                    help="(do you) want to include a pheatmap? TRUE or FALSE")
parser$add_argument("-d", "--groupsincoldendrogram", type="character", 
                    help="Number of groups in the col Dendrogram from the pheatmap")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="output folder where you want to store your results")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

#############################
## Reading or preparing the inputs
#############################
mymatrix <-read.table( file=args$matrix, stringsAsFactors = FALSE , check.names = FALSE)
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toys_merged_quantile_norm_by_batch.txt", stringsAsFactors = FALSE, check.names = FALSE)
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Merged/Exp_Mat_GSE115989_MNHK.tsv", stringsAsFactors = FALSE, check.names = FALSE)
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Merged/Exp_Mat_MK_GSE113342LE_GSE115989RJ_MajaL_GSE89880.txt", stringsAsFactors = FALSE, check.names = FALSE)
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swipped_LipCon2-3_for_NorChl2-3.txt", stringsAsFactors = FALSE, check.names = FALSE) 

annotdf <-read.table( file=args$annotation, stringsAsFactors = FALSE , header=TRUE )
# annotdf <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toys_merged_annotations.tsv", stringsAsFactors = FALSE )
# annotdf <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Merged/GSE115989_MNHK_AnnotFile.tsv", stringsAsFactors = FALSE )
# annotdf <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Merged/Annot_MK_GSE113342_GSE115989_ML_GSE89880.tsv", stringsAsFactors = FALSE , header=TRUE )
# annotdf <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/annotation_lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_09_Rformat.tsv", stringsAsFactors = FALSE , header=TRUE)

code_path <- args$code
# code_path <- "/media/rmejia/mountme88/code/Processing_an_RNAexpression_matrix_from_hyb_n_scanning/"

label <- args$label # label <- "your_title" # label <- "lipidosis_RNA_16_STAR_fC_edgeR_matrix_Log_transformed_swipped_LipCon23-NorChl23"
your_main_groups <- args$maingroups # your_main_groups <- "group"

myperplexitynumber <- args$perplexity
mode(myperplexitynumber) <- "numeric"
# myperplexitynumber <- 3

PCA_point_size <- args$pcapointsize
mode( PCA_point_size ) <- "numeric"
# PCA_point_size <- 3

Include_pheatmap_TRUEorFALSE <- args$whatapheatmat
# Include_pheatmap_TRUEorFALSE <- "TRUE"

Number_of_groups_in_the_cols_dendrogram_from_the_heatmap <- args$groupsincoldendrogram
# Number_of_groups_in_the_cols_dendrogram_from_the_heatmap <- 5

outputfolder <- args$outputfolder
#  outputfolder <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/toy_merged_quantiles_normalized_per_batch/"
#  outputfolder <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/GSE115989_MNHK/"
#  outputfolder <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/Annot_MK_GSE113342_GSE115989_ML_GSE89880"
#  outputfolder <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_idk_if_they_are_shuffled/swipped_LipCon23-NorChl23_log_transformed_Heatmap"

dir.create(outputfolder, recursive = TRUE)
outputfolder <- normalizePath(outputfolder)

code_path <- normalizePath(code_path)

##############################
## The program starts
#############################
if(all(colnames( mymatrix) ==  annotdf$Unique_ID) == TRUE ){ 
  print("your annotation and colnames match")
}
if(all(colnames( mymatrix) ==  annotdf$Unique_ID) != TRUE ){ 
  print("ERROR: your annotation and colnames DON´T match")
  break()
}

####################
### Filtering the matix, eliminating the zero rows
####################
RowPositions_with_all_zeros <- rowSums(mymatrix) == 0 # row positions with only zeros
matrix_excluding_AllZeroRows <- mymatrix[!RowPositions_with_all_zeros, ] # eliminating those rows

#####################
## Annotation for the plots
####################
annot_4_plotting_pca <- annotdf 
annot_4_plotting_pca[ , your_main_groups ] <- as.factor( annot_4_plotting_pca[ , your_main_groups ] )

source( paste0( code_path ,"/libraries/" , "matrix_N_annotdf_2_melteddf.R") )
melted_matrix_excluding_AllZeroRows <- matrix_N_annotdf_2_melteddf( matrix_excluding_AllZeroRows  , annot_4_plotting_pca )
# melted_expmat <- matrix_N_annotdf_2_melteddf( mymatrix , annot_4_plotting_pca )

############
## graphs
############
source(paste0( code_path,"/libraries/","PCA-choosing-poing-size_box_density_pheatmap_tsnes_plots.R") )
PCA_box_density_pheatmap_tsnes_plots( paste0(  outputfolder,"/PCA_2D" )  ,
                                      matrix_excluding_AllZeroRows ,
                                      annot_4_plotting_pca ,
                                      melted_matrix_excluding_AllZeroRows ,
                                      paste0( label ) ,
                                      myperplexitynumber ,
                                      your_main_groups, PCA_point_size , 
                                      Number_of_groups_in_the_cols_dendrogram_from_the_heatmap,
                                      Include_pheatmap_TRUEorFALSE
)

# write.table(expmat_log2, file= paste0(path_to_your_matrix ,"_log2.tsv") , sep="\t", row.names = TRUE, col.names = TRUE )