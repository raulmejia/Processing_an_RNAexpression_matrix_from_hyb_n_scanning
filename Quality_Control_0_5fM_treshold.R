# This script do the Quality control for Ncounter data
# It calculates 2 standard deviations of the probe POS_E (That is supposed to be 0.5fM https://www.nanostring.com/wp-content/uploads/2020/12/Gene_Expression_Data_Analysis_Guidelines.pdf)
# and adds that quantitiy to the mean of the Negative probes, 2*sd(POS_E)+mean(Negative_Controls) = the cut off
# The structure of your matrix should be 
#           HNGC(or other id) sample1 sample2 ...
# rowname1  0.42    45
# rowname2  1       0
# rowname3  45      0
# NEG_Prob1 0       12
# POS_E     2       1

## Notes
## Make it flexible in order to change the name of the probe (POS_E) or other and the number of SD 
##
## Should I substact the treshould value to the final matrix or not?
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
your_feature <- "POS_E"
inputdf <-read.table( file=args$inputfile, stringsAsFactors = FALSE )
# inputdf <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toy_for_treshold.txt", stringsAsFactors = FALSE)

path2save <- args$outputfile
#  path2save <- "/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/toy_for_treshold.txt_fM5.5_treshold.tsv"

# getting the mean from the negative controls
mean_each_NEG_probe <- apply(inputdf[grep("NEG_",rownames(inputdf)),], 1, mean)
mean_across_NEG_probes <- mean(mean_each_NEG_probe)
treshold <- mean_across_NEG_probes + 2*sd(inputdf[your_feature,])

logical_matrix <- inputdf > treshold
rows_with_more_than_2_samples_above_the_treshold <- apply(logical_matrix,1,sum) > 1

prunned_matrix <- inputdf[ rows_with_more_than_2_samples_above_the_treshold , ]

###########################
#   Saving the results  ###
###########################
write.table( prunned_matrix, file = path2save , row.names = TRUE, sep="\t", col.names = TRUE )

###########################
#   Saving the lost rows ###
###########################
lost_rows <- setdiff( rownames(inputdf) , rownames(prunned_matrix) )
write.table( lost_rows, file = paste0(path2save,"_lost_rows") , row.names = TRUE, sep="\t", col.names = TRUE )


