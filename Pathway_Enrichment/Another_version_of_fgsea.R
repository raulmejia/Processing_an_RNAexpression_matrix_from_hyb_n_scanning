# GSE try me
# https://bioinformatics-core-shared-training.github.io/RNAseq_May_2020_remote/html/06_Gene_set_testing.html
library(tidyverse)
library(fgsea)
#library("org.Hs.eg.db")


df_with_my_ranking_path <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/DESeq2_RNAseq_lipidosis_vs_Normal_First_sheet_gseaGudenas_format_untilColumn8_Nonovel-X-pseudogene-d-uncharacterized_nopoints.tsv"
df_with_my_ranking <- read.table(file=df_with_my_ranking_path, header = TRUE,sep = "\t")
S4table <- read.table(file = df_with_my_ranking_path, sep='\t', header=TRUE) %>%
  filter(Gene.Symbol != "")

gene_list = S4table[ , Column_with_your_rankins ]
names(gene_list) = S4table$Gene.Symbol
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]



rankData <- gene_list


load("/media/rmejia/mountme88/Common_and_Virgin_Data/Common_Data_Across_Projects/Pathways/GO/human_c5_v5p2.rdata")
ls()


pathwaysH <- Hs.c5

fgseaRes <- fgsea(pathwaysH, 
                 rankData, 
                 minSize=15, 
                 maxSize = 500, 
                 nperm=1000)

fgseaRes %>% 
  arrange(desc(abs(NES))) %>% 
  top_n(10, -padj)
