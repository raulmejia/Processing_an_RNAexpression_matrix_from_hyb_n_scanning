############
## this script receives information about gene expressión and pathway definitions and retrieves
##   the consensus regulated pathways (gage and fgsea) according to the function of Brian Gudenas  https://bioinformaticsbreakdown.com/how-to-gsea/
##  Inputs:
##    A)  A data frame (Note that it should be in csv format)
##          It should contain your rankins
##          And a column "Gene.Symbol" to give the name of the genes to this ranking
##    B)  The name of the column with your rankings = Column_with_your_rankings
##    C)  Path to Code
##    D)  

## your column of interest should be named like "Gene.Symbol"
## the column of your rankings should be named "DESeq2.Log2.Fold.Change"

# lets look at https://bioinformatics-core-shared-training.github.io/RNAseq_May_2020_remote/html/06_Gene_set_testing.html

if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("dplyr")) {
  BiocManager::install("dplyr", ask =FALSE)
  library("dplyr")
}
if (!require("fgsea")) {
  BiocManager::install("fgsea", ask =FALSE)
  library("fgsea")
}
if (!require("gage")) {
  BiocManager::install("gage", ask =FALSE)
  library("gage")
}
if (!require("data.table")) {
  BiocManager::install("data.table", ask =FALSE)
  library("data.table")
}
if (!require("ggplot2")) {
  BiocManager::install("ggplot2", ask =FALSE)
  library("ggplot2")
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
parser$add_argument("-d", "--dataframe", type="character", 
                    help="path to your data frame in csv format")
parser$add_argument("-r", "--columnrankings", type="character", 
                    help="name of the column that contains your rankings")
parser$add_argument("-c", "--code", type="character", 
                    help="path to your code")
parser$add_argument("-l", "--label", type="character", 
                    help="label (small one) for your results")
parser$add_argument("-p", "--pathwaysfile", type="character", 
                    help="Pathway file (in .gmt) format ")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="output folder where you want to store your results")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

##############################
### Reading the input data
##############################

code_path <- args$code
# code_path <- "/media/rmejia/mountme88/code/Processing_an_RNAexpression_matrix_from_hyb_n_scanning/Pathway_Enrichment/"
code_path <- normalizePath(code_path)

source(paste0(code_path,"/GSEAgudenas.R"))
source(paste0(code_path,"/Fgsea_npermSimple_intersection_gage.R"))

df_with_my_ranking_path <- args$dataframe
# df_with_my_ranking_path <- "/home/rmejia/Downloads/toy/pone.0145322.s007_4.csv"
##### df_with_my_ranking_path <- "/home/rmejia/Downloads/toy/pone.0145322.s007.tsv"
##### df_with_my_ranking_path <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/DESeq2_RNAseq_lipidosis_vs_Normal_First_sheet_gseaGudenas_format_untilColumn8_Nonovel-X-pseudogene-d-uncharacterized_nopoints.tsv"
##### df_with_my_ranking_path <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/DESeq2_Output_ALL_Deleted_Mar-Sep.tsv"
##### df_with_my_ranking_path <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/DESeq2_Output_ALL_Mar_ToMACH-Sep_ToSEPT-Dec_ToDEC.tsv"

#####df_with_my_ranking <- read.table(file=df_with_my_ranking_path, sep = "\t")
df_with_my_ranking <- read.table(file=df_with_my_ranking_path, header = TRUE,sep = "\t")

GO_filepath <- args$pathwaysfile
# GO_filepath <- "/media/rmejia/mountme88/Common_and_Virgin_Data/Common_Data_Across_Projects/Pathways/MsigDB/C5_GeneOntology_Biological_Process/c5.go.bp.v7.4.symbols.gmt"

Column_with_your_rankins <- args$columnrankings
# Column_with_your_rankins <- "DESeq2.Log2.Fold.Change"

mysmalllabel4title<- args$label
# mysmalllabel4title <- "WT--vs--DT"

output_folder <- args$outputfolder
# output_folder <- "/home/rmejia/Downloads/toy/Results/DIL/ALL_DESeq_Mar-Sep-Dec-Rescued"
dir.create( output_folder, recursive = TRUE )
output_folder <- normalizePath(output_folder)

mypval <- 0.05
###############################
### The program starts
###############################
# S4table <- read.csv(file = df_with_my_ranking_path, skip=1) %>%  filter(Gene.Symbol != "")
 S4table <- read.table(file = df_with_my_ranking_path, sep='\t', header=TRUE) %>%
  filter(Gene.Symbol != "")

gene_list = S4table[ , Column_with_your_rankins ]
names(gene_list) = S4table$Gene.Symbol
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

#  myGO = fgsea::gmtPathways(GO_filepath)
# 
# fgRes <- fgsea::fgsea(pathways = myGO, 
#                       stats = gene_list,
#                       minSize=15,
#                       maxSize=600,
#                       nperm=10000)


fgRes <- fgsea::fgsea(pathways = myGO, 
                      stats = gene_list,
                      minSize=15, # Do it lower than 10 ? Check the ABCA4 Path
                      maxSize=600) %>% as.data.frame() %>%  dplyr::filter(padj < !!mypval)

fgRes <- fgsea::fgsea(pathways = myGO, 
                      stats = gene_list,
                      minSize=15, # Do it lower than 10 ? Check the ABCA4 Path
                      nPermSimple = 10000,
                      maxSize=600) %>% as.data.frame() %>%  dplyr::filter(padj < !!mypval)
fgRes_sorted_by_padj <- fgRes[order(fgRes$padj, decreasing = TRUE) ,]

dim(fgRes_sorted_by_padj)

#Rescue the MARCH and SEPT

# fgRes_Multilevel <- fgseaMultilevel(
#   pathways = myGO, 
#   stats = gene_list,
#   minSize=15, # Do it lower than 10 ? Check the ABCA4 Path
#   maxSize=600
#   ) %>%  dplyr::filter(padj < !!mypval)

gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,600))
str(gaRes)
attr(gaRes$greater,"dimnames")
head(gaRes$less, n= 20)
gaRes$stats

res = Fgsea_npermSimple_intersection_gage( gene_list, GO_filepath , pval = 0.05 , mysmalllabel4title )  
res = GSEAgudenas( gene_list, GO_filepath , pval = 0.05 , mysmalllabel4title )
# res = GSEAgudenas( gene_list[c(1:100,18227:18327)], GO_filepath , pval = 0.05 , mysmalllabel4title )
dim(res$Results)



?fgsea

pdf( file = paste0(output_folder ,"/" , mysmalllabel4title ,".pdf") ,
     width = 10, height = 7)
print(res$Plot)
dev.off()

# Matrix_to_save <- res$Results[ , c("pathway","pval","padj","ES","Enrichment")]
Matrix_to_save <- res$Results
save_matrix_path <- paste0(output_folder ,"/" , mysmalllabel4title ,"_matrix.tsv" )
write.table( as.matrix(Matrix_to_save) , file=save_matrix_path , sep="\t", row.names = FALSE )


save_rds_path <- paste0(output_folder ,"/" , mysmalllabel4title ,".RDS" )
saveRDS( res, save_rds_path)
