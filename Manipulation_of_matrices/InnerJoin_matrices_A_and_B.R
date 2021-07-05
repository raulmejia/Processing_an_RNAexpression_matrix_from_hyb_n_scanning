# This script merge two matrices (Inner Join according SQL operations https://es.wikipedia.org/wiki/Sentencia_JOIN_en_SQL)
# They should have the same ids in the rows (in order to make exact match)
#          sample1 sample2 ...
# rowname1  0.42    45
# rowname2  1       0
# rowname3  45      0
# NEG_Prob1 0       12
# POS_E     2       1
#
# Example of use: Rscript /path/Merge_2_matrices_by_intersected_rows.R -A /path/path2/some_expmatrix.tsv -B /path1/path3/another_matrix2.tsv -o /path3/merged_matrices.tsv
## Notes:
## Make it flexible in order to change the thresholds (number of samples) and minimum value 
############################## 
## Required libraries
##############################
if (!require("argparse")) {
  BiocManager::install("argparse", dependencies = TRUE)
  library("argparse")
}
if (!require("dplyr")) {
  BiocManager::install("dplyr", dependencies = TRUE)
  library("dplyr")
}
if (!require("ggplot2")) {
  BiocManager::install("ggplot2", dependencies = TRUE)
  library("ggplot2")
}
if (!require("limma")) {
  BiocManager::install("limma", dependencies = TRUE)
  library("limma")
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
parser$add_argument("-A", "--infileA", type="character", 
                    help="input file A")
parser$add_argument("-B", "--infileB", type="character", 
                    help="input file B")
parser$add_argument("-o", "--outputfile", type="character", 
                    help="output file where you want to store your results")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
# print some progress messages to stderr if "quietly" wasn't requested

#############################
## The program starts
#############################
inputA <-read.table( file=args$infileA, stringsAsFactors = FALSE , check.names=FALSE)
# inputA <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toy_for_treshold.txt", stringsAsFactors = FALSE)

inputB <-read.table( file=args$infileB, stringsAsFactors = FALSE , check.names=FALSE)
# inputB <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toy_for_treshold_B.txt", stringsAsFactors = FALSE)

path2save <- args$outputfile
#  path2save <- "/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toys_merged.txt"

common_rows  <- intersect( rownames(inputA ) , rownames(inputB) )
posA <- which( rownames(inputA ) %in% common_rows )
A_filtered <- inputA[ posA, ]
posB <- which( rownames(inputB ) %in% common_rows )
B_filtered <- inputB[ posB, ]
result_df <- cbind( A_filtered , B_filtered )
  

###########################
#   Saving the results  ###
###########################
write.table( result_df , file = path2save , row.names = TRUE, sep="\t", col.names = TRUE )

######################################
#   Some proportions of the result ###
######################################
print("percentage of rows from A included in the merged matrix")
print( dim(result_df)[1] / dim(inputA)[1] )

print("percentage of rows from B included in the merged matrix")
print( dim(result_df)[1] / dim(inputB)[1] )
