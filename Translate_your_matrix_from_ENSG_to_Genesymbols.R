library("gprofiler2")

gensymbols <-read.table("/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29-switched-to_work_colsymb-n-gene-deleted_original_order.tsv_log2.tsv", sep="\t", header =TRUE)

head(gensymbols)

rownames(gensymbols)
queryENSG  <- gsub( "\\.[:graph:]*","" , rownames(gensymbols))

gprofiler2::gconvert(query = queryENSG[1:3], organism = "hsapiens",
         target="HGNC", mthreshold = Inf, filter_na = TRUE)









gprofiler2::gconvert(query = rownames(gensymbols)[1:3], organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)


converted <- gconvert(rownames(gensymbols) , organism = "hsapiens", target = "HGNC",
         region_query = F, numeric_ns = "", mthreshold = Inf,
         filter_na = T, df = T)


converted <- gconvert(rownames(gensymbols) , organism = "hsapiens", target = "ENSG",
                      region_query = F, numeric_ns = "", mthreshold = Inf,
                      filter_na = T, df = T)


?gProfileR::gconvert()


if (!require("biomaRt")) {
  BiocManager::install("biomaRt", dependencies = TRUE)
  library("biomaRt")
}


ensembl = useMart("ensembl")
ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

hgnc_symbol_list <- rownames(gensymbols)
list_with_queries <- list()
for(k in 1:length(hgnc_symbol_list) ){
  biomresult <- biomaRt::getBM(attributes='ensembl_gene_id', 
                               filters = 'hgnc_symbol', 
                               values = hgnc_symbol_list[k], 
                               mart = ensembl)
  if(length(biomresult[,1]) == 0){
    list_with_queries[[k]] <- "NA"
  }
  if(length(biomresult[,1]) == 1){
    list_with_queries[[k]] <- biomresult[1,1]
  }
  
}

retrievable <- data.frame(hgnc_symbol_list , as.character(list_with_queries))
colnames(retrievable) = c("gene_symbols", "ensembl_ids")

gconvert(query = c("REAC:R-HSA-3928664", "rs17396340", "NLRP1"), organism = "hsapiens",
         target="ENSG", mthreshold = Inf, filter_na = TRUE)


gconvert(query = c("REAC:R-HSA-3928664", "rs17396340", "NLRP1"), organism = "hsapiens",
         target="ENSG", mthreshold = Inf, filter_na = TRUE)





