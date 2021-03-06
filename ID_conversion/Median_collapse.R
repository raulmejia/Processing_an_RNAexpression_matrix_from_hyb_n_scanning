# This script takes a expression matrix and performs a median collapsing of your repeated IDS 
#
# The structure of your matrix should be 
#           HNGC(or other id) sample1 sample2 ...
# rowname1  HNGCA             0.42    45
# roename2  HNGCB             1       0
# rowname3  HNGCB             45      0

# Example of use:
# Rscript Median_collapse.R \
#  -i /Path/to/input_matrix.tsv \
#  -o /Path/to/output_matrix.tsv

## Notes
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
parser$add_argument("-i", "--inputfile", type="character", 
                    help="input file with your gene list in genesymbols")
parser$add_argument("-o", "--outputfile", type="character", 
                    help="output file where you want to store your results")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
# print some progress messages to stderr if "quietly" wasn't requested

#############################
## The program starts
#############################
#############################
## 1) Reading the data
#############################
inputdf <-read.table( file=args$inputfile, stringsAsFactors = FALSE )
# inputdf <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_HNGC_added.tsv", stringsAsFactors = FALSE)

path2save <- args$outputfile
#  path2save <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_HNGC_added_median_collapsed.tsv"

#############################
## 2) Checking if it pertinent to proceed 
#############################
if(sum(duplicated(inputdf[,1])) == 0 ){
  print( "your matrix doesn't have replicates so it is already collapsed, there is nothing to do in median collapse" )
  write.table(inputdf, path2save, sep = "\t", row.names = TRUE, col.names = TRUE)
  break()
}

#############################
## 3) Separating the IDs and the numeric matrix
#############################
your_ids_to_collapse <- inputdf[,1]
the_matrix<- inputdf[,-1]

#############################
## 4) Identifying the positions of the duplicated and non-duplicated IDs
#############################
duplicated_ids <- inputdf[,1][duplicated(inputdf[,1])]
duplicated_positions <- which(inputdf[,1] %in% duplicated_ids) # Identifying the positions where the duplicates lay

positions_non_duplicated <- setdiff( 1:dim(inputdf)[1] , duplicated_positions) # Identifying the complementary positions ( where the duplicates haven't lain)
df_of_non_duplicated <- inputdf[ positions_non_duplicated , ] # Data frame with NO duplicates

df_of_duplicated <- inputdf[ duplicated_positions , ] # Data frame of the duplicates
# sort(table(df_of_duplicated[,1]  )) # Do you want to see a table with your duplicates? 

#############################
## 4) Preparing the data and the function to obtain the median per ID across the samples
#############################
df_of_duplicated[,1] <- as.factor(df_of_duplicated[,1]) # Factor to split the data frames
df_of_duplicated_HNGC_splitted <- split(df_of_duplicated,df_of_duplicated[,1])  # Getting a list of "mutually exclusive dataframes". Each data frame contains all the repetitions of the "non-unique" ID 

median_col_by_data_frame <-function(somedf ){ # function to get the median by cols from a dataframe
  df_for_median <- somedf[-1]
  median_vector<-apply(df_for_median,2,median)
  #median_vector_plus_id <- c(unique(somedf[1,]) , median_vector)
  return( median_vector )
}
# median_col_by_data_frame(df_of_duplicated_HNGC_splitted[[1]])

collapsed_by_median_lists <- lapply( df_of_duplicated_HNGC_splitted , median_col_by_data_frame ) # Applying the function to the list of "mutually exclusive dataframes"
df_of_duplicated_now_collapsed_by_median <- do.call(rbind, collapsed_by_median_lists) # Transforming the list of DFs to a single Matrix

#############################
## 5) Joining the Data Frame with unique IDs and the DataFrame of the duplicated IDs (Now collapsed by median)
#############################
rownames( df_of_non_duplicated ) <- df_of_non_duplicated[ , 1 ] # We can assign this rownames bc there is no duplicated
df_of_non_duplicated_no_id_colum <- df_of_non_duplicated[ , -1] # Now that the IDs can be used as rownames we can deleted the first column

df_median_collapsed <- rbind(df_of_non_duplicated_no_id_colum,  df_of_duplicated_now_collapsed_by_median) # binding the two DFs
df_median_collapsed <- df_median_collapsed[ order(rownames( df_median_collapsed ) ) , ] # ordering the rownames

###########################
# 6)  Saving the results  ###
###########################
write.table( df_median_collapsed, file = path2save , row.names = TRUE, sep="\t" )
