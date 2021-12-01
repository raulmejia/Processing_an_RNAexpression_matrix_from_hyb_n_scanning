# Internal gender through X and Y expression levels

# FilteringMatrix has quotes and it is separated by tabs
# Example

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

filteringmatrix <-read.table( file=args$matrixtofilter , stringsAsFactors = FALSE , check.names = FALSE)
# filteringmatrix <- read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Testing_swapping_rearraegments/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swapped_LipCon2_NorChl2_Shifted2Left_LipChl3_LipCon3_NorChl3.txtAnnotations_Extracting_X_n_Y_and_annotating.tsv", stringsAsFactors = FALSE , header=TRUE)

namecol_fromMF <- args$namecolfrommatrixfilter
# namecol_fromMF <-"ensembl_gene_id"

annotdf <-read.table( file=args$annotation, stringsAsFactors = FALSE , header=TRUE )
# annotdf <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Data/Annotations/annotation_lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_09_Rformat.tsv", stringsAsFactors = FALSE , header=TRUE)

your_main_groups <-"Unique_ID"

code_path <- args$code
# code_path <- "/media/rmejia/mountme88/code/Processing_an_RNAexpression_matrix_from_hyb_n_scanning/"

label <- args$label # label <- "your_title" #
#your_main_groups <- args$maingroups # your_main_groups <- "group"

outputfolder <- args$outputfolder
#  outputfolder <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Testing_swapping_rearraegments"

dir.create(outputfolder, recursive = TRUE)
outputfolder <- normalizePath(outputfolder)
code_path <- normalizePath(code_path)

##############################
## The program starts
#############################
# if(all(colnames( mymatrix) ==  annotdf$Unique_ID) == TRUE ){ 
#   print("your annotation and colnames match")
# }

# sum(rownames(mymatrix) %in% filteringmatrix[,namecol_fromMF] )

rownames_mymatrix<-  rownames(mymatrix)
Common_Positions <- which(str_extract(rownames_mymatrix, regex("[:alnum:]*")) %in% filteringmatrix[,namecol_fromMF])
rowstouse <- rownames_mymatrix[Common_Positions]

matrix_prunned <-  mymatrix[ which(rownames(mymatrix) %in% rowstouse ) , ]

# if(all(colnames( mymatrix) ==  annotdf$Unique_ID) != TRUE ){ 
#   print("ERROR: your annotation and colnames DON´T match")
#   break()
# }

####################
### Filtering the matix, eliminating the zero rows
####################
RowPositions_with_all_zeros <- rowSums(matrix_prunned) == 0 # row positions with only zeros
matrix_excluding_AllZeroRows <- matrix_prunned[!RowPositions_with_all_zeros, ] # eliminating those rows

#####################
## Annotation for the plots
####################
annot_4_plotting_pca <- annotdf 
annot_4_plotting_pca[ , your_main_groups ] <- as.factor( annot_4_plotting_pca[ , your_main_groups ] )

expmat_log2 <- log2( matrix_excluding_AllZeroRows +1 )

source( paste0( code_path ,"/libraries/" , "matrix_N_annotdf_2_melteddf.R") )
melted_expmat_log2 <- matrix_N_annotdf_2_melteddf(expmat_log2 , annot_4_plotting_pca )





