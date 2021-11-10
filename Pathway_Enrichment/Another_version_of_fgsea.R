# GSE try me
# https://bioinformatics-core-shared-training.github.io/RNAseq_May_2020_remote/html/06_Gene_set_testing.html
library(tidyverse)
library(fgsea)
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("clusterProfiler")) {
  BiocManager::install("clusterProfiler", ask =FALSE)
  library("clusterProfiler")
}

BiocManager::install("clusterProfiler")
if (!require("org.Hs.eg.db")) {
  BiocManager::install("org.Hs.eg.db", ask =FALSE)
  library("org.Hs.eg.db")
}
library(clusterProfiler)
#library("org.Hs.eg.db")

library(clusterProfiler)
##############################
### Info given by the user
##############################
Path2yourDataFrame <-"/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/DESeq2_Output_ALL.txt"
Pathway_data_base_GMTformat_path <- "/media/rmejia/mountme88/Common_and_Virgin_Data/Common_Data_Across_Projects/Pathways/MsigDB/C5_GeneOntology_Biological_Process/c5.go.bp.v7.4.symbols.gmt"
Column_with_your_rankins <- "DESeq2.Log2.Fold.Change"
ColName_of_your_genes <- "Gene.Symbol"

###########################
### The program starts here
###########################
yourDataFrame <-read.table(file = Path2yourDataFrame , sep="\t" , header=TRUE) # Reading the matrix 

# yourDataFrame_No_blank_spaces <- yourDataFrame %>% filter( Gene.Symbol != "" )
# dim(yourDataFrame_No_blank_spaces)
# yourDataFrame_No_blank_spaces[ 3069:3072, ]

# Prefiltering
yourDataFrame_No_blank_spaces_neitherMarOrSep <-
  yourDataFrame %>% filter( Gene.Symbol != "" ) %>% filter(!grepl("Mar ", Gene.Symbol)) %>% filter(!grepl("Sep ", Gene.Symbol))
# Deleting empty rows and errors induced by Exel SEPT12 -> sep12

FilteredDF <- yourDataFrame_No_blank_spaces_neitherMarOrSep # changing to a shorter name

myranked_list = FilteredDF[ , Column_with_your_rankins ] # Getting ready the ranked list
names(myranked_list) = FilteredDF[ ,ColName_of_your_genes ]
myranked_list = sort( myranked_list , decreasing = TRUE)
myranked_list = myranked_list[ !duplicated( names( myranked_list ) ) ] 

res = GSEAgudenas( myranked_list, GO_filepath , pval = 0.05 , mysmalllabel4title )

dim(res$Results)
res$Results
print(res$Plot)

########################################################


load("/media/rmejia/mountme88/Common_and_Virgin_Data/Common_Data_Across_Projects/Pathways/GO/human_c5_v5p2.rdata")
load("/media/rmejia/mountme88/Common_and_Virgin_Data/Common_Data_Across_Projects/Pathways/GO/mouse_c5_v5p2.rdata")
ls()

"Mm.c5"
pathwaysH <- Hs.c5

fgseaRes <- fgsea(pathwaysH, 
                 rankData, 
                 minSize=15, 
                 maxSize = 500, 
                 nperm=1000)

fgseaRes <- fgsea(Mm.c5, 
                  rankData, 
                  minSize=15, 
                  maxSize = 500, 
                  nperm=1000)



fgseaRes %>% 
  arrange(desc(abs(NES))) %>% 
  top_n(10, -padj)

library(clusterProfiler)
library(enrichplot)


# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# reading in data from deseq2
df = read.csv("expression_set_full-gsea.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)
write.csv(df, file ="dfff.csv")


# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")





















