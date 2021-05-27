# This script tries to correct the batch effect throug SVA 
# 
# Annotation file:
#   One column of your annotation file should be named "Unique_ID" and must match in exaclty the same order with the rownames of the expression matrix
#   Another column describe the batches you can choose the name of such a column
#
# Structure of the expression matrix:
#          sample1 sample2 ...
# rowname1  0.42    45
# rowname2  1       0
# rowname3  45      0
# NEG_Prob1 0       12
# POS_E     2       1
#
# Put a general description of the input files in the README of the repository
###################################
# The files who contain the pretended matrices should use "." Decimal insted of ","
# The annotation file should contain columns with exactly the following names: "Morph_cat_Andre","Histology_number","Scan_ID","Biopsy_year","Morphological_Categories"
###################################
#### 0) loading and/or installing required libraries
###################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("sva")) {
  BiocManager::install("sva", ask =FALSE)
  library("sva")
}
if (!require("limma")) {
  BiocManager::install("limma", ask =FALSE)
  library("limma")
}
if (!require("reshape2")) {
  install.packages("reshape2", ask =FALSE)
  library("reshape2")
}
if (!require("preprocessCore")) {
  BiocManager::install("preprocessCore", ask =FALSE)
  library("preprocessCore")
}
if (!require("affy")) {
  BiocManager::install("affy", ask =FALSE)
  library("affy")
}

############################## 
## Data given by the user
##############################
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]" )
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-m", "--matrix", type="character", 
                    help="input matrix")
parser$add_argument( "-c", "--code", type="character" , 
                     help="path to your code file" )
parser$add_argument("-a", "--annotation", type="character", 
                    help="input annotation file")
parser$add_argument("-b", "--batches", type="character", 
                    help="Column describing your batches in your annotation file")
parser$add_argument("-o", "--outputfile", type="character", 
                    help="output file where you want to store your results")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
# print some progress messages to stderr if "quietly" wasn't requested

#############################
## The program starts
#############################
inmatrix <-read.table( file=args$matrix, stringsAsFactors = FALSE )
# inmatrix <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toys_merged_quantile_norm_by_batch.txt", stringsAsFactors = FALSE)

annot <-read.table( file=args$annotation, stringsAsFactors = FALSE )
# annot <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toys_merged_annotations.tsv", stringsAsFactors = FALSE)

path2save <- args$outputfile
#  path2save <- "/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toys_merged_quantile_norm_by_batch_Batch_correction_through_SVA.txt"

code_path <- args$code
# code_path <- "/media/rmejia/mountme88/code/Processing_an_RNAexpression_matrix_from_hyb_n_scanning/"
code_path <- normalizePath(code_path)

batches_col <- args$batches
# batches_col <- "group"

#############################
## The program starts
#############################
#####################
### Check the inputs
#####################
if(all(colnames(inmatrix) ==  annot$Unique_ID) == TRUE ){ 
  print("your annotation and colnames match")
}
if(all(colnames(inmatrix) ==  annot$Unique_ID) != TRUE ){ 
  print("ERROR: your annotation and colnames DON´T match")
  break()
}

############
#### running SVA
###########
mypheno <- as.data.frame(annot[,batches_col])
rownames(mypheno) <- annot[,"Unique_ID"]
colnames(mypheno) <- c("batch")

mypheno$batch <- as.factor(mypheno$batch)
batch <- mypheno$batch
modcombat <- model.matrix(~1, data=mypheno)

combat_qnormBsep = ComBat(dat = inmatrix , batch=batch , mod=modcombat, par.prior=TRUE, prior.plots=FALSE)


###########################
#   Saving the results  ###
###########################
write.table( combat_qnormBsep  , file = path2save , row.names = TRUE, sep="\t", col.names = TRUE )
