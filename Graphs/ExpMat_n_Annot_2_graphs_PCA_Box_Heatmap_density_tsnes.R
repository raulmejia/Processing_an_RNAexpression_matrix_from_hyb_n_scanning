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
## Notes

# Please note that the columns in the annotation file should no contain numbers otherwise the tsnes complain
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
#  mymatrix <-read.table(file="/data/...Exp_Mat_GSE113342_Tubulus.tsv", stringsAsFactors = FALSE, check.names = FALSE)

annotdf <-read.table( file=args$annotation, stringsAsFactors = FALSE, header = TRUE )
# annotdf <-read.table(file="/data/.../Annotation_GSE113342_Only_TubuloInterstitium.tsv", stringsAsFactors = FALSE, , header = TRUE)

code_path <- args$code
# code_path <- "/data/.../Ncounter_RCC_processing/"

label <- args$label # label <- "As_downloaded_tubulointerstitium"
your_main_groups <- args$maingroups # your_main_groups <- "Tissue"

myperplexitynumber <- args$perplexity
mode(myperplexitynumber) <- "numeric"
# myperplexitynumber <- 3

PCA_point_size <- args$pcapointsize
mode( PCA_point_size ) <- "numeric"
#PCA_point_size <- 3

outputfolder <- args$outputfolder
#  outputfolder <- "/data/.../As_Downloaded/"

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

## Plotting pcas
annot_4_plotting_pca <- annotdf
annot_4_plotting_pca[ , your_main_groups ] <- as.factor( annot_4_plotting_pca[ , your_main_groups ] )

source( paste0( code_path ,"/libraries/" , "matrix_N_annotdf_2_melteddf.R") )
meltedrawdata <- matrix_N_annotdf_2_melteddf( mymatrix , annotdf )

############
## graphs
############
source(paste0( code_path,"/libraries/","PCA-choosing-poing-size_box_density_heatmap_tsnes_plots.R") )
PCA_box_density_heatmap_tsnes_plots( paste0(  outputfolder,"/PCA_2D" )  ,
                             mymatrix ,  annot_4_plotting_pca , meltedrawdata , paste0( label ) , myperplexitynumber , your_main_groups, PCA_point_size )
