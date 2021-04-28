
exp_matrix <- Raw_expmat
annotdf <- annot_4_plotting_pca
melteddf <- meltedrawdata
  label4title <- data_label
  result_dir <- path_Results_directory
  
PCA_box_density_plots <- function( result_dir, exp_matrix, annotdf, melteddf, label4title  ){
  exp_matrix_T <- t(exp_matrix)
  dir.create( result_dir, recursive = TRUE )
  pdf( file=paste0(result_dir,"/" , label4title ,".pdf"),
       width = 10, height = 7)
  print(autoplot( prcomp( exp_matrix_T ), data = annotdf, colour= 'group') +
          ggtitle(paste(label4title )))
  q <- ggplot( melteddf , aes(Unique_ID, value, fill=group))
  print( q + geom_boxplot( )+ ggtitle(paste( label4title )) )
  
  
  plot1 <- ggplot(melteddf, aes(x=value, fill = group, y = Unique_ID)) + 
    geom_density_ridges() + 
    ggtitle(label4title)
  print(plot1)
  #pheatmap( exp_matrix ,cutree_cols = 5,  annotation_col  = annotdf[,c("Histology_number","Morphological_Categories","Scan_ID" )], fontsize = 4) 
  #pheatmap( exp_matrix ,cutree_cols = 5, col = brewer.pal(  length(table(annotdf[,"Morphological_Categories"])) ,"Set3"),  annotation_col  = annotdf[,c("group" )], fontsize = 4) 
  
  tsne_model_1 = Rtsne(  t(exp_matrix)  , check_duplicates=FALSE, pca=TRUE, perplexity=5, theta=0.5, dims=2)
  d_tsne_1 = as.data.frame( tsne_model_1$Y )
  
  d_tsne_1_simplecols <- cbind(d_tsne_1, annotdf$group)
  mytsneplot_Nocolors <- ggplot(d_tsne_1_simplecols, aes(x=V1, y=V2 )) +
    geom_point(size=2.5) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle( paste("t-SNE",label4title) ) +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) +
    scale_colour_brewer(palette = "Set3")
  print(mytsneplot_Nocolors)

  dev.off()  
}