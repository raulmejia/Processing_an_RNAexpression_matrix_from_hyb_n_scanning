# This script converts the ENSBL to genesymbols from a matrix 
# The final numbers after the point are the version
# https://m.ensembl.org/Help/Faq?id=488

## Notes
## this program erase the numbers after the point in ENS codes while left the .XXX_PAR_region
## untouched if you want to erase it you can substitute for this chunk of code: ensemblids_noversion <- sub( '\\.[0-9]*$', '', ensemblids_version )
## for this one: ensemblids_noversion <- sub( '\\.[[:graph:]]*$', '', ensemblids_version ) # Deletes PAR label but creates replicated 
## Nevertheless an undesirable effect of this chance is the duplicity in ENS codes before the point

############################## 
## Required libraries
##############################
if (!require("gprofiler2")) {
  BiocManager::install("gprofiler2", dependencies = TRUE)
  library("gprofiler2")
}
if (!require("argparse")) {
  BiocManager::install("argparse", dependencies = TRUE)
  library("argparse")
}
if (!require("robustbase")) {
  BiocManager::install("robustbase", dependencies = TRUE)
  library("robustbase")
}
if (!require("dplyr")) {
  BiocManager::install("dplyr", dependencies = TRUE)
  library("dplyr")
}
if (!require("parallel")) {
  BiocManager::install("parallel", dependencies = TRUE)
  library("parallel")
}
if (!require("MASS")) {
  BiocManager::install("MASS", dependencies = TRUE)
  library("MASS")
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
inputmatrix <-read.table( file=args$inputfile, stringsAsFactors = FALSE )
# inputmatrix <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29-switched-to_work_colsymb-n-gene-deleted_original_order.tsv_log2.tsv", stringsAsFactors = FALSE)
# inputmatrix <-read.table(file="/data/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29-switched-to_work_colsymb-n-gene-deleted_original_order.tsv_log2.tsv", stringsAsFactors = FALSE)

path2save <- args$outputfile
#  path2save <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_23_HNGC_collapsed_by_median.tsv"
#  path2save <- "/data/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_23_HNGC_collapsed_by_median.tsv"

print("printing your arguments")
print(args$inputfile)
print(args$outputfile)

ensemblids_version <- rownames( inputmatrix )
ensemblids_noversion <- sub( '\\.[0-9]*$', '', ensemblids_version )

# getting the map between ESN (ensembl) and HGNC
table_HGNCids <- gprofiler2::gconvert( ensemblids_noversion , organism = "hsapiens", target="HGNC" )

# Checking the table with the converted Ids
# Identifying the input ids that mapped to more than 1 HNGC
duplicatedinputids <- table_HGNCids$input[ duplicated( table_HGNCids$input ) ] # The duplicated ENSG
duplicatedinputids_positions <- which(table_HGNCids$input %in% duplicatedinputids) # 
inputids_that_map_to_more_than_1_HGNC <- table_HGNCids[ duplicatedinputids_positions , ] # ENSG that map to more than 1 HGNC in a table  
# save this table besides the results
write.table( inputids_that_map_to_more_than_1_HGNC , file = paste0(path2save,"_ENSG_that_map_to_more_than_1_HNGC.txt") , row.names = TRUE, sep="\t" )

# Identifying the HGNCs that received more than one map
duplicatedHGNCsids <- table_HGNCids$name[ duplicated( table_HGNCids$name ) ]
duplicatedHGNCs_positions <- which(table_HGNCids$name %in% duplicatedHGNCsids)
table_of_duplicated_HGNCs <- table_HGNCids[duplicatedHGNCs_positions, ]

#######
## Bulding the matrix with the annotations (HGNC)
#######
input_ids_and_converted_ids <- table_HGNCids[ , c(2,4)]

# creating the matrix to fill out
exp_values_part <- matrix( rep(0, dim( input_ids_and_converted_ids)[1] * dim(inputmatrix)[2] ) , ncol = dim(inputmatrix)[2] ) 
converted_ids_uncollapsed_table <- cbind( input_ids_and_converted_ids, exp_values_part) 
colnames(converted_ids_uncollapsed_table)[ 3:dim(converted_ids_uncollapsed_table)[2] ] <- colnames(inputmatrix  )

##########
## filling the extended dataframe (now the df has more more rows bc of the repeated ENSG that mapped to more than 1 HNGC)
##########
  give_me_the_values_from_the_mastrix <- function( mycharacter){
    position_ENSG_in_inputmat_develop <- grep(  mycharacter, rownames(inputmatrix) )
    return(inputmatrix[ position_ENSG_in_inputmat_develop, ])
  }

  start_time_mcapply <- Sys.time()  # Taking the time to the parallel use of the defined function
  numCores <- detectCores()
  expressed_values_per_input_row <- mclapply(converted_ids_uncollapsed_table[,"input"], give_me_the_values_from_the_mastrix , mc.cores = numCores -1) # mclapply is the star for this parallelization I used all the cores but 1
  end_time_mcapply <- Sys.time()
  print(end_time_mcapply - start_time_mcapply)

  df_numeric_expression_with_results_gconvert <- do.call(rbind,expressed_values_per_input_row  ) # Binding the results from the previous code
  
  HGNC_df <- input_ids_and_converted_ids
  if( all(sub( '\\.[0-9]*$', '', rownames( df_numeric_expression_with_results_gconvert ) ) == HGNC_df$input) == TRUE ) {
    print("Successful final check, the rownames of the translated numeric matrix and the result from gprofiler2::gconvert matches")
  }
  if( all(sub( '\\.[0-9]*$', '', rownames( df_numeric_expression_with_results_gconvert ) ) == HGNC_df$input) == FALSE ) {
    print("ERROR: final check, the rownames of the translated numeric matrix and the result from gprofiler2::gconvert doesn't match")
  }
  
  rownames(HGNC_df) <- rownames( df_numeric_expression_with_results_gconvert) # preparing the HNGC df to paste to the final df
  HGNC_df <- HGNC_df[-grep("input",colnames(HGNC_df))]
  colnames(HGNC_df)[grep("target",colnames(HGNC_df))] <- "HNGC"

  result_data_frame <- cbind(HGNC_df, df_numeric_expression_with_results_gconvert )# adding the HNGC to the final df
  
###################
# Saving the results  
###################
  write.table(file=path2save,result_data_frame, sep="\t", row.names = TRUE, col.names = TRUE  )
