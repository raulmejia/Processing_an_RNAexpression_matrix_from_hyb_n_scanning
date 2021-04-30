# Input a numeric matrix, an annotation file, and a row of the matrix that we want to generate the box plot



# Or take only the differential expressed genes for the clustering

# Delete generation 2?
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

