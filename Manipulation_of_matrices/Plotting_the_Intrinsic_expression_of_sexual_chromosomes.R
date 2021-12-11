# Internal gender through X and Y expression levels

# FilteringMatrix has quotes and it is separated by tabs
# Example
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", ask =FALSE)
  library("RColorBrewer")
}
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("dplyr")) {
  BiocManager::install("dplyr", ask =FALSE)
  library("dplyr")
}
if (!require("ggplot2")) {
  BiocManager::install("ggplot2", ask =FALSE)
  library("ggplot2")
}
if (!require("reshape2")) {
  install.packages("reshape2", ask =FALSE)
  library("reshape2")
}
if (!require("stringr")) {
  install.packages("stringr", ask =FALSE)
  library("stringr")
}
if (!require("pheatmap")) {
  BiocManager::install("pheatmap", ask =FALSE)
  library("pheatmap")
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
parser$add_argument("-f", "--matrixtofilter", type="character", 
                    help="path to your filtering matrix")
parser$add_argument("-n", "--namecolfrommatrixfilter", type="character", 
                    help="what's the col name to use from your filtering matix")
parser$add_argument("-c", "--code", type="character", 
                    help="path to your code")
parser$add_argument("-l", "--label", type="character", 
                    help="label to your results")
parser$add_argument("-g", "--maingroups", type="character", 
                    help="the name of your column to correct / make intrabatch normalization")
parser$add_argument("-p", "--parregions", type="character", 
                    help="Include Pseudoautosomal regions or not TRUE/FALSE")
parser$add_argument("-s", "--specialsample", type="character", 
                    help="The name of a special sample from which you want an individual boxplot")
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
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Data/Rearrangements_of_the_Source/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swapped_LipCon2_NorChl2_Shifted2Left_LipChl3_LipCon3_NorChl3.txt", stringsAsFactors = FALSE, check.names = FALSE) 
#  mymatrix <-read.table(file="/data/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Data/Rearrangements_of_the_Source/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swapped_LipCon2_NorChl2_Shifted2Left_LipChl3_LipCon3_NorChl3.txt", stringsAsFactors = FALSE, check.names = FALSE) 
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Data/Source/lipidosis_RNA_16_STAR_fC_edgeR_matrix.txt", stringsAsFactors = FALSE, check.names = FALSE) 

filteringmatrix <-read.table( file=args$matrixtofilter , stringsAsFactors = FALSE , check.names = FALSE)
# filteringmatrix <- read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Testing_swapping_rearraegments/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swapped_LipCon2_NorChl2_Shifted2Left_LipChl3_LipCon3_NorChl3.txtExtracting_Chr_X_and_annotating.tsv", stringsAsFactors = FALSE , header=TRUE)
# filteringmatrix <- read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Testing_swapping_rearraegments/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swapped_LipCon2_NorChl2_Shifted2Left_LipChl3_LipCon3_NorChl3.txtExtracting_Chr_Y_and_annotating.tsv", stringsAsFactors = FALSE , header=TRUE)
# filteringmatrix <- read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Testing_swapping_rearraegments/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swapped_LipCon2_NorChl2_Shifted2Left_LipChl3_LipCon3_NorChl3.txtExtracting_Chr_Y_and_annotating.tsv", stringsAsFactors = FALSE , header=TRUE)

namecol_fromMF <- args$namecolfrommatrixfilter
# namecol_fromMF <-"ensembl_gene_id"

annotdf <-read.table( file=args$annotation, stringsAsFactors = FALSE , header=TRUE )
# annotdf <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Data/Annotations/annotation_lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_09_Rformat.tsv", stringsAsFactors = FALSE , header=TRUE)
# annotdf <-read.table(file="/data/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Data/Annotations/annotation_lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_09_Rformat.tsv", stringsAsFactors = FALSE , header=TRUE)

#### !# Do I need this ???  your_main_groups <-"Unique_ID"

code_path <- args$code
# code_path <- "/media/rmejia/mountme88/code/Processing_an_RNAexpression_matrix_from_hyb_n_scanning/"
# code_path <- "/data/code/Processing_an_RNAexpression_matrix_from_hyb_n_scanning/"

label <- args$label # label <- "your_title" #
# your_main_groups <- args$maingroups # your_main_groups <- "group"

PAR_regions <- args$parregions
# PAR_regions <- FALSE

special_sample <- args$
# special_sample <- "LipCon2"

outputfolder <- args$outputfolder
#  outputfolder <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Testing_swapping_rearraegments"
#  outputfolder <- "/data/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Testing_swapping_rearraegments"

dir.create(outputfolder, recursive = TRUE)
outputfolder <- normalizePath(outputfolder)
code_path <- normalizePath(code_path)

##############################
## The program starts
#############################

####################
### Applying the PAR criteria
####################
filteringmatrix_afterPAR <- filteringmatrix # to compare after the application of the filter
if( PAR_regions == FALSE){
  rows_according_PAR_criteria <- is.na(filteringmatrix$PAR1) & is.na( filteringmatrix$PAR2 )
  filteringmatrix_afterPAR <- filteringmatrix[rows_according_PAR_criteria,] 
}

print("dim of the filtering matrix before and after PAR criteria")
dim(filteringmatrix) 
dim(filteringmatrix_afterPAR)

####################
### Applying the Filtering matrix
####################
#rownames_mymatrix <-  rownames(filteringmatrix_afterPAR ) # Selecting the rows out of your Matrix according to the given filtered matrix
Common_Positions <- which( str_extract( rownames(mymatrix), regex("[:alnum:]*") ) %in% filteringmatrix_afterPAR [,namecol_fromMF] ) # Getting the positions for selecting the rows out of your matrix
#rowstouse <- rownames_mymatrix[Common_Positions] # 
matrix_prunned <-  mymatrix[ Common_Positions  , ] 
#matrix_prunned <-  mymatrix[ which(rownames( filteringmatrix_afterPAR  ) %in% rowstouse ) , ] 

dim( matrix_prunned)[1] ; dim( filteringmatrix)[1] # Proportion between the prunned matrix and the filtering matrix
print ("dimension of the filtering matrix given as input and after applying the criteria")

####################
### Eliminating the Rows with only zeros
####################
RowPositions_with_all_zeros <- rowSums(matrix_prunned) == 0 # row positions with only zeros
matrix_excluding_AllZeroRows <- matrix_prunned[!RowPositions_with_all_zeros, ] # eliminating those rows

print ("dimension of the matrix after excluding all zero rows") ; dim(matrix_excluding_AllZeroRows) 

#####################
## Annotation for the plots
####################
annot_4_plotting_pca <- annotdf 
annot_4_plotting_pca[ , your_main_groups ] <- as.factor( annot_4_plotting_pca[ , your_main_groups ] )

expmat_log2 <- log2( matrix_excluding_AllZeroRows +1 )

source( paste0( code_path ,"/libraries/" , "matrix_N_annotdf_2_melteddf.R") )
melted_expmat_log2 <- matrix_N_annotdf_2_melteddf(expmat_log2 , annot_4_plotting_pca )

melted_expmat <- matrix_N_annotdf_2_melteddf(matrix_excluding_AllZeroRows , annot_4_plotting_pca )


pdf(file= paste0(outputfolder,"/", label) , height = 7, width = 7 )

g <- ggplot(melted_expmat_log2, aes( x= Unique_ID, y=value) ) +
  geom_boxplot()
g  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#g+ coord_flip()  + geom_jitter(shape=16, position=position_jitter(0.2))

gem <- ggplot(melted_expmat , aes( x= Unique_ID, y=value) ) +
  geom_boxplot()
gem  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


####
# Expecial sample to visualize

melted_expmat_log2_special_sample <- melted_expmat_log2 %>% filter( Unique_ID == special_sample  )
f <- ggplot(melted_expmat_log2_special_sample, aes( x=  variable, y=value) ) +
  geom_boxplot()
f  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

#RowPositionsmymatrix_with_all_zeros<- rowSums(mymatrix) == 0 # row positions with only zeros
#Mymatrix_excluding_AllZeroRows <- mymatrix[!RowPositionsmymatrix_with_all_zeros, ] # eliminating those rows

#log2Mymatrix_excluding_AllZeroRows <- log2( Mymatrix_excluding_AllZeroRows +1 )
#heatmap( as.matrix( log2Mymatrix_excluding_AllZeroRows ))
#cal_z_score <- function(x){ # Defining the funtion for the z-transformation
#  (x - mean(x)) / sd(x)
#}
#expmat_zscore <- t( apply( log2Mymatrix_excluding_AllZeroRows , 1, cal_z_score ) )

#pheatmap( expmat_zscore )
