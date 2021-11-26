
# https://datavizpyr.com/heatmaps-in-r-with-pheatmap-package/
# https://r-charts.com/correlation/pheatmap/
# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("pheatmap")) {
  BiocManager::install("pheatmap", ask =FALSE)
  library("pheatmap")
}
if (!require("dendextend")) {
  BiocManager::install("dendextend", ask =FALSE)
  library("dendextend")
}

mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_swipped_LipCon2-3_for_NorChl2-3.txt", stringsAsFactors = FALSE, check.names = FALSE)


# Delete the all zero rows
RowPositions_with_all_zeros <- rowSums(mymatrix) == 0
matrix_excluding_AllZeroRows <- mymatrix[!RowPositions_with_all_zeros, ]

log2mat <- log2( matrix_excluding_AllZeroRows +1 )
dim(mymatrix)
dim(log2mat)

boxplot(t(log2mat))

pheatmap( mymatrix[1:10000,] )
pheatmap( log2mat[1:10000,] )
pheatmap( log2mat ) 

cal_z_score <- function(x){ # Defining the funtion for the z-transformation
  (x - mean(x)) / sd(x)
}

log2mat_zscore <- t( apply( log2mat, 1, cal_z_score ) )
pheatmap( log2mat_zscore[1:10000,])

my_hclust_gene <- hclust( dist ( log2mat_zscore ) , method = "ward.D2" )
# my_hclust_gene <- hclust( dist ( log2mat_zscore ) , method = "complete" )

my_hclust_sample <- hclust( dist ( t(log2mat_zscore[1:100,]) ) , method = "ward.D2" )
as.dendrogram(my_hclust_sample) %>%
  plot(vertical = TRUE)

my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)