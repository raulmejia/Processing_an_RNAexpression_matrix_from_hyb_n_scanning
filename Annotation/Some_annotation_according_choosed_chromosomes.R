###################################
#### This script receives a matrix with ENSMBL ids as rownames and samples as columns it
####    should retreieve a dataframe with the ensembls as rownames and annotations (chromosome_name","start_position","end_position","hgnc_symbol) 
####    for the X and Y chromosomes it will say if it is a Pseudoautosomal region or not) it will also retrieve a DF with the Pseduoautosomal region information
####    Input:
#LipChl1  LipChl2 ContUnt1  ContUnt2  LipUnt1
#ENSG00000215601 0.42 45  25  45  65
#ENSG00000215603 1 0 2 3 1
#ENSG00000231341 45 0 12  23  2 
#ENSG00000231436 0 2 1 12 15
#ENSG00000234795 2 5 6 4 8
####      
####     Output:
####         A list inlcuding this dataframe, a data fram with the PAR positions and two lists with the PAR genes
#ensembl_gene_id chromosome_name start_position end_position hgnc_symbol
#ENSG00000167393  X 333933  386955  PPP2R3B
#ENSG00000225661  X 1008503 1010101 RPL14P5
####      
#### Example:
####     
####Rscript /Path/to/Script/Some_annotation_according_choosed_chromosomes.R \
#### -m /Path/to/Expression_Matrix/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swapped_LipCon2_NorChl2_Shifted2Left_LipChl3_LipCon3_NorChl3.txt \
#### -c /Path/2/your/chromosomes/of/interest/chromXnY.txt \
#### -l Some_label_for_the_titles_and_file_names_of_your_results \
#### -o /Path/for/the/output/  
####
####  To do: 
####     Save in an R object the PAR genes and the PAR table  
####     Do another program that preserve the genes and just give annotations
###################################
#### 0) loading and/or installing required libraries
################################### 
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("org.Hs.eg.db")) {
  BiocManager::install("org.Hs.eg.db", ask =FALSE)
  library("org.Hs.eg.db")
}
if (!require("biomaRt")) {
  BiocManager::install("biomaRt", ask =FALSE)
  library("biomaRt")
}
if (!require("tidyverse")) {
  install.packages("tidyverse", ask =FALSE)
  library("tidyverse")
}
if (!require("httr")) {
  install.packages("httr", ask =FALSE)
  library("httr")
}

library(httr)
set_config(config(ssl_verifypeer = 0L))
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
parser$add_argument("-c", "--chomosomes", type="character", 
                    help="path to your expression matrix")
parser$add_argument("-l", "--label", type="character", 
                    help="label to your results")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="output folder where you want to store your results")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

#############################
## Reading or preparing the inputs
#############################

# mymatrixpath <-"/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Data/Rearrangements_of_the_Source/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swapped_LipCon2_NorChl2_Shifted2Left_LipChl3_LipCon3_NorChl3.txt"
mymatrixpath <- args$matrix
mymatrix <-read.table( file = mymatrixpath , stringsAsFactors = FALSE , check.names = FALSE)

label <- args$label #  label <- "Annotations_Extracting_X_n_Y_and_annotating.tsv"

# Choosen_Chr_path <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Testing_swapping_rearraegments/chromXnY.txt"
Choosen_Chr_path <- args$chomosomes
Choosen_Chr <- read.table( file=Choosen_Chr_path, stringsAsFactors = FALSE , check.names = FALSE, header = FALSE )
# myChosedChromosomes <- c("X","Y")
Choosen_Chr <- Choosen_Chr$V1

outputfolder <- args$outputfolder
# outputfolder <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Testing_swapping_rearraegments/"
dir.create( outputfolder, recursive = TRUE )
outputfolder <- normalizePath( outputfolder )
#dirname( ) #basename( )

##############################
## The program starts
#############################
##############################
## Preparing the Mart
#############################
ensembl <-                                ## fully specified mart
  useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#head(listMarts(), 3)                      ## list the marts
#head(listDatasets(useMart("ensembl")), 3) ## mart datasets
#head(listFilters(ensembl), 20)             ## filters
myFilter <- "chromosome_name"
myValues <- Choosen_Chr
myAttributes <- c("ensembl_gene_id","chromosome_name","start_position","end_position","hgnc_symbol" )
res <- getBM(attributes =  myAttributes, filters =  myFilter,
             values =  myValues, mart = ensembl)

##############################
## adding the PAR annotation
#############################
PAR1 <- sort( c( "PLCXD1","GTPBP6","PPP2R3B","SHOX","CRLF2","CSF2RA","IL3RA","SLC25A6"
                 ,"ASMTL","P2RY8","AKAP17A","ASMT","DHRSX","ZBED1","CD99","XG" ))
PAR2 <- c( "IL9R", "SPRY3","VAMP7" )

res$PAR1 <- rep( NA , dim( res )[ 1 ] )
PAR1positions <- which(res$hgnc_symbol %in% PAR1)
res$PAR1[PAR1positions] <- "PAR1"

res$PAR2 <- rep( NA , dim( res )[ 1 ] )
PAR2positions <- which(res$hgnc_symbol %in% PAR2)
res$PAR2[PAR2positions] <- "PAR2"

#saving the DF
final_path_to_save <-paste0( outputfolder,"/", basename( mymatrixpath ) ,label )
write.table( res, file= final_path_to_save , sep="\t", col.names = TRUE, row.names = TRUE, quote=TRUE )

##############################
## The dataframe with the PAR annotation
#############################
PositionsPARs <- data.frame( version=c( "GRCh38","GRCh38","GRCh38","GRCh38","GRCh37","GRCh37","GRCh37","GRCh37" ) ,
                             chromosome_name = c( "X","X","Y","Y","X","X","Y","Y" ) ,
                             start_position= c(10001,155701383,10001,56887903,60001,154931044,10001,59034050  ) ,
                             end_position = c(2781479 ,156030895,2781479,57217415,2699520,155260560,2649520,59363566 ),
                             chromosome_name = c( "PAR1","PAR2","PAR1","PAR2","PAR1","PAR2","PAR1","PAR2" ) )

write.table( PositionsPARs, file = paste0( outputfolder, "/" , basename(mymatrixpath ) , label , "PARS_info.tsv" ) )
