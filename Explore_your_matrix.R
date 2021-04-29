###################################
#### Author: Raúl Mejía
#### The aim of this program is to explore your expression matrix, taking a glimpse of the normalization & batch effects  
###################################
# one column of your annotation file should have the name "group" and other "Unique_ID"
# the column names of your expression matrix should match with the rows entries of the "group" column in your annotation df (in thar exact order)
# This script calls the function PCA_box_density_plots_group_Treatment_Cell_line that plot the group Treatment and Cell_line columns from your annotation file (those columns should be in your annotation file)
###################################
#### 0) loading and/or installing required libraries
###################################
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
if (!require("sva")) {
  BiocManager::install("sva", ask =FALSE)
  library("sva")
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
if (!require("cowplot")) {
  install.packages("cowplot", ask =FALSE)
  library("cowplot")
}
if (!require("preprocessCore")) {
  BiocManager::install("preprocessCore", ask =FALSE)
  library("preprocessCore")
}
if (!require("affy")) {
  BiocManager::install("affy", ask =FALSE)
  library("affy")
}
if (!require("oligo")) {
  BiocManager::install("oligo", ask =FALSE)
  library("oligo")
}
if (!require("Rtsne")) {
  BiocManager::install("Rtsne", ask =FALSE)
  library("Rtsne")
}
#if (!require("M3C")) {
#  BiocManager::install("M3C", ask =FALSE)
#  library("M3C")
#}
#if (!require("tidyverse")) {
#  BiocManager::install("tidyverse", ask =FALSE)
#  library("tidyverse")
#}

###################################
#### Data given by the user
###################################
# myargs <- commandArgs(trailingOnly = TRUE)
path_to_your_matrix <- myargs[1]
# path_to_your_matrix <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29_to_work_colsymb-n-gene-deleted.tsv"
# path_to_your_matrix <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29-switched-to_work_colsymb-n-gene-deleted_original_order.tsv"

path_to_your_annotation_file <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/annotation_lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_09_Rformat.tsv"
# path_to_your_annotation_file <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/annotation_lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_09.tsv"

Code_path <- "/media/rmejia/mountme88/code/Processing_an_RNAexpression_matrix_from_hyb_n_scanning"  
# Path where are the rest of your scripts

path_Results_directory <-"/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/Results_from_the_exploratory_analysis"

# data_label<- "lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29_to_work_colsymb-n-gene-deleted"
# data_label<- "switched_org_orderlipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_29_to_work_colsymb-n-gene-deleted"

colname_4_intra_batch_normalization <- "group" # the name of your column to correct
# please don't use spaces or parenthesis/brackets in the names of your columns

###################################
#### Normalize your paths
###################################
Code_path <- normalizePath(Code_path)
path_Results_directory  <- normalizePath( path_Results_directory  )

###################################
#### Creating your result folders if they doesn't exist yet
###################################
dir.create(path_Results_directory , recursive = TRUE)

###################################
#### Reading the annotation table and the table that cointains the expression data
###################################
Raw_expmat <- read.table( path_to_your_matrix  , sep = "\t", header = TRUE)
annot <- read.table( path_to_your_annotation_file , sep = "\t", header = TRUE )

#####################
# Annotation object for plotting pcas
####################
annot_4_plotting_pca <- annot
annot_4_plotting_pca[ , "group" ] <- as.factor( annot_4_plotting_pca[ , "group" ] )
colnames(Raw_expmat) ==  annot_4_plotting_pca$Unique_ID

# loading the function to melt (reshape) the data to preparation for ggplot2 functions
source( paste0( Code_path,"/libraries/","matrix_N_annotdf_2_melteddf.R") )
meltedrawdata <- matrix_N_annotdf_2_melteddf( Raw_expmat , annot )

########################################
########
########    0. Exploratory
########
########################################
#########
### Visualize your the Raw data
########
source(paste0( Code_path,"/libraries/","PCA_box_density_plots.R") )
PCA_box_density_plots_group_Treatment_Cell_line(  paste0( path_Results_directory,"/Exploratory" )  ,
                        Raw_expmat ,  annot_4_plotting_pca , meltedrawdata , paste0( data_label, "Data_as_given" ))

########################################
########
########    1. Preprocessing
########
########################################
##################
## log 2 transformation
##################
expmat_log2 <- log2( Raw_expmat +1)
melted_expmat_log2 <- matrix_N_annotdf_2_melteddf(expmat_log2 , annot )

# Visualize the data log2 transformed data
PCA_box_density_plots_group_Treatment_Cell_line(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        expmat_log2 ,  annot_4_plotting_pca , melted_expmat_log2 , paste0( data_label, "_Log2" ))






###################################
#### Normalization intra batch
###################################
source(paste0( Code_path,"/","Matrix_2_list_of_sub_matrices.R") )
list_of_submatrices <- Matrix_2_list_of_sub_matrices( table$Scan_ID , expmat_log2 )

################################
###### quantile normalization  (normalizeQuantiles) batch Separated
################################
list_splitd_qnorm <- lapply( list_of_submatrices  ,  normalizeQuantiles)
mat_qnorm_sep_by_batch <- do.call(cbind, list_splitd_qnorm)

# All in once q norm
qnorm_all_of_once <- normalizeQuantiles( expmat_log2 )

#### visualize your normalized data
melted_mat_qnorm_sep_by_batch <- matrix_N_annotdf_2_melteddf( mat_qnorm_sep_by_batch , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        mat_qnorm_sep_by_batch ,  annot_4_plotting_pca , melted_mat_qnorm_sep_by_batch , paste( data_label, "Log2_Qnorm_intra" ))

melted_qnorm_all_of_once <- matrix_N_annotdf_2_melteddf( qnorm_all_of_once , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        qnorm_all_of_once ,  annot_4_plotting_pca , melted_qnorm_all_of_once  , paste( data_label, "Log2_Qnorm_all_inonce" ))


###########
### cyclic loess normalization 
############
################################
###### cyclic loess normalization batch Separated
################################
list_splitd_CLnorm <- lapply( list_of_submatrices  ,  function(x) {normalizeCyclicLoess(x, method = "affy")} )
mat_cyclicloess_norm_sep_by_batch <- do.call(cbind, list_splitd_CLnorm)

# all in once
#normalizeCyclicLoess( expmat_log2)
expmat_log2_cyclic_loess_AO <- normalizeCyclicLoess( expmat_log2, method = "fast")

# visualization
melted_mat_cyclicloess_norm_sep_by_batch  <- matrix_N_annotdf_2_melteddf( mat_cyclicloess_norm_sep_by_batch  , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        mat_cyclicloess_norm_sep_by_batch  ,  annot_4_plotting_pca , melted_mat_cyclicloess_norm_sep_by_batch , paste( data_label, "Log2_CLnorm_intra" ))

melted_expmat_log2_cyclic_loess_AO <- matrix_N_annotdf_2_melteddf( expmat_log2_cyclic_loess_AO  , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        expmat_log2_cyclic_loess_AO  ,  annot_4_plotting_pca , melted_expmat_log2_cyclic_loess_AO , paste( data_label, "Log2_CLnorm_all_inonce" ))


#############
## Looking for batch effect
#############
## Pheno object for combat
#pheno <- annot[, c("Morph_cat_Andre","Scan_ID")]
pheno <- annot[, c("Morphological_Categories","Scan_ID")]
colnames(pheno) <- c("subgroups","batch")
rownames(pheno) <- annot[,"Unique_ID"] 

# Batch effect with combat
pheno$batch <- as.factor(pheno$batch)
batch <- pheno$batch
modcombat <- model.matrix(~1, data=pheno)
combat_qnormBsep = ComBat(dat = mat_qnorm_sep_by_batch , batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
combat_qnormAllofOnce = ComBat(dat = qnorm_all_of_once , batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
combat_cyclic_loess_BatchSep = ComBat(dat =  mat_cyclicloess_norm_sep_by_batch , batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
combat_cyclic_loess_AO = ComBat(dat = expmat_log2_cyclic_loess_AO , batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#modcombat<-model.matrix(~subgroups, data=pheno)
#combat_mydata<-ComBat(dat= t(mymatrix), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

### Visualization
meltedcombat_qsepareted <- matrix_N_annotdf_2_melteddf( combat_qnormBsep , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        combat_qnormBsep ,  annot_4_plotting_pca , meltedcombat_qsepareted , paste( data_label, "Log2_QnormBatchSep_combat" ))

meltedcombat_qnormAllinOne <- matrix_N_annotdf_2_melteddf( combat_qnormAllofOnce , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        combat_qnormAllofOnce ,  annot_4_plotting_pca , meltedcombat_qnormAllinOne , paste( data_label, "Log2_QnormInBatch_combatAcross" ))

melted_combat_cyclic_loess_BatchSep <- matrix_N_annotdf_2_melteddf( combat_cyclic_loess_BatchSep , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        combat_cyclic_loess_BatchSep ,  annot_4_plotting_pca , melted_combat_cyclic_loess_BatchSep , paste( data_label, "Log2_CycLoessBatchSep_combat" ))

meltedcombat_cyclic_loess_AO <- matrix_N_annotdf_2_melteddf( combat_cyclic_loess_AO , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        combat_cyclic_loess_AO ,  annot_4_plotting_pca , meltedcombat_cyclic_loess_AO , paste( data_label, "Log2_CycLoessAO_combat" ))

write.table( combat_qnormBsep, paste0( path_Results_directory,"/Preprocessing","/",paste0( data_label, "Log2_QnormBatchSep_combat" ),".tsv" ), sep="\t" )
write.table( combat_cyclic_loess_BatchSep, paste0( path_Results_directory,"/Preprocessing","/",paste0( data_label, "Log2_CycLoessBatchSep_combat" ),".tsv" ), sep="\t" )

## MA plot Checking qualiy controls 
limma::plotMA(  expmat_log2, array=3, main = c("MA plot exp Mat" , "array_number",3)  )
limma::plotMA( combat_qnormAllofOnce, array=3)
limma::plotMA( combat_qnormAllofOnce, array=1)
limma::plotMA(  expmat_log2, array=1)


# After DEAnalysis
## https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/plotMA

# RNA degradation with "AffyRNAdeg"


mydesign <- model.matrix( ~ 0 + annot$Morphological_Categories)
colnames(mydesign) <- gsub("annot\\$Morphological_Categories","",colnames(mydesign))
contrast.matrix <- makeContrasts(paste0("Normal","-",setdiff(colnames(mydesign), "Normal")[1]),
                                 paste0("Normal","-",setdiff(colnames(mydesign), "Normal")[2]),
                                 paste0("Normal","-",setdiff(colnames(mydesign), "Normal")[3]),
                                 paste0("Normal","-",setdiff(colnames(mydesign), "Normal")[4]),
                                 paste0("Normal","-",setdiff(colnames(mydesign), "Normal")[5]),levels=mydesign)
fit <- lmFit( combat_qnormBsep , mydesign )
fita1c <-contrasts.fit(fit , contrast.matrix)
fita1c<-eBayes(fita1c)



fit<-lmFit(a1c[,1:200],design)
contrast.matrix<-makeContrasts(Normal-Medium, Normal-High, Medium-High, levels=design)
fita1c<-contrasts.fit(fit,contrast.matrix)
fita1c<-eBayes(fit2)
fita1c<-eBayes(fita1c)



table( annot$Morphological_Categories )


dim(combat_cyclic_loess_AO)


combat_qnormBsep

# Create design matrix for leukemia study
testchr<- as.character(annot$Morphological_Categories)
testchr[ testchr == "Normal"  ] <- "Normal"
testchr[ testchr != "Normal"  ] <- "disease"
test2df   <- data.frame(testchr)
colnames(test2df) <- "Mycat"

design <- model.matrix(~Mycat, data = test2df)

# Count the number of samples modeled by each coefficient
colSums(design)

# Load package
library(limma)

# Fit the model
fit <- lmFit(combat_qnormBsep, design)

# Calculate the t-statistics
fit <- eBayes(fit)

# Summarize results
results <- decideTests(fit[, "Mycatnormal"])
results <- decideTests(fit[, "MycatNormal"])

#options(digits=2) 

genes<- topTable(fit, coef=2, n=60, adjust="BH") 
summary(results)


#---------------------------------
Morph_cat_chr <- as.character(annot$Morphological_Categories)
log_Norm_necr_w_cell <- (grepl("Normal",Morph_cat_chr ) | grepl("necrosis_w_cell_cresc",Morph_cat_chr ) )
Normal_nec_w_cresc_df <- data.frame(Morph_cat_chr [log_Norm_necr_w_cell] )
colnames(Normal_nec_w_cresc_df) <- "Mycat"

design_crec_w <- model.matrix(~Mycat, data = Normal_nec_w_cresc_df)

# Fit the model
fit_crec_w <- lmFit( combat_qnormBsep[,log_Norm_necr_w_cell], design_crec_w )

# Calculate the t-statistics
fit_crec_w <- eBayes(fit_crec_w)

# Summarize results
results <- decideTests(fit_crec_w[, "MycatNormal"])

#options(digits=2) 

genes <- topTable(fit, coef=2, n=60, adjust="BH") 
genes[,c("logFC","adj.P.Val")]
summary(results)


#--------------------------------
Morph_cat_chr <- as.character(annot$Morphological_Categories)
Normal_nec_w_cresc[ Normal_nec_w_cresc == "Normal"  ] <- "Normal"
log_fibrocell_Crescent <- (grepl("Normal",Morph_cat_chr ) | grepl("fibrocell_Crescent",Morph_cat_chr ) )
Normal_fibrocell_Crescent_df <- data.frame(Morph_cat_chr[log_fibrocell_Crescent] )
colnames(Normal_fibrocell_Crescent_df) <- "Mycat"

design_ibrocell_Crescent <- model.matrix(~Mycat, data = Normal_fibrocell_Crescent_df)

# Fit the model
fit_fibrocell_Crescent <- lmFit( combat_qnormBsep[,log_fibrocell_Crescent] , design_ibrocell_Crescent )

# Calculate the t-statistics
fit_fibrocell_Crescent <- eBayes(fit_fibrocell_Crescent)

# Summarize results
results_fibrocell_Crescent <- decideTests(fit_fibrocell_Crescent[, "MycatNormal"])

#options(digits=2) 

genes_fibrocell_Crescent <- topTable(fit_fibrocell_Crescent , coef=2, n=60, adjust="BH") 
genes_fibrocell_Crescent[,c("logFC","adj.P.Val")]
summary(results)


#-----------------------------

Morph_cat_chr <- as.character(annot$Morphological_Categories)
log_necrosis_only <- (grepl("Normal",Morph_cat_chr ) | grepl("necrosis_only",Morph_cat_chr ) )
Normal_necrosis_only_df <- data.frame(Morph_cat_chr[log_necrosis_only] )
colnames(Normal_necrosis_only_df) <- "Mycat"

design_necrosis_only <- model.matrix(~Mycat, data = Normal_necrosis_only_df)

# Fit the model
fit_necrosis_only <- lmFit( combat_qnormBsep[,log_necrosis_only] , design_necrosis_only )

# Calculate the t-statistics
fit_necrosis_only <- eBayes(fit_necrosis_only)

# Summarize results
results_necrosis_only <- decideTests(fit_necrosis_only[, "MycatNormal"])

#options(digits=2) 

genes_necrosis_only <- topTable(fit_necrosis_only , coef=2, n=60, adjust="BH") 

genes_necrosis_only[,c("logFC","adj.P.Val")]
summary( results_necrosis_only )






