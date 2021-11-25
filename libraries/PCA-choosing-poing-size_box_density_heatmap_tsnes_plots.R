#####################
# this program plots a PCA, boxplot, denstity plot and Tsne for each of the meaningful variables* that are included in your annoation file
# your annotation file needs a column "Unique_ID"
#
# Note: 
# Search what should be the apropriate perplexity number for the tsnes
####################

# exp_matrix <-   mymatrix
# annotdf <- annot_4_plotting_pca
# melteddf <-meltedrawdata
# label4title <- paste0( label )
# result_dir <- paste0( outputfolder,"/PCA_2D" )
# myperplexity <- 3 # (default was 5)
# group_4_tsnes <- your_main_groups
# PCA_point_size <-3

# Please note that the columns in the annotation file should no contain numbers otherwise the tsnes complain

#############
##A required function
#############
#  someannotdf <- annotdf 
#  colname_from_your_annordf <- group_4_tsnes 
# k =1
colors_4_plotDensities <- function( someannotdf , colname_from_your_annordf ){
  group_uniquenames_length <- length( unique( someannotdf[,colname_from_your_annordf] ) ) # length of the uniquenames
  rainbow_colors <- rainbow( group_uniquenames_length ) # color to substitute those names
  rainbow_colors_4_plot <- rep( "NA" , length( someannotdf[,colname_from_your_annordf] )   )
  for( k in 1:group_uniquenames_length) {
    positions_K <- someannotdf[,colname_from_your_annordf] == unique( someannotdf[,colname_from_your_annordf] )[k]
    rainbow_colors_4_plot[positions_K] <- rainbow_colors[k] 
  }
  return(rainbow_colors_4_plot)
}

####
# The program starts
####
PCA_box_density_heatmap_tsnes_plots <- function( result_dir, exp_matrix, annotdf, melteddf, label4title, myperplexity, group_4_tsnes, PCA_point_size){
  exp_matrix_T <- t(exp_matrix)
  dir.create( result_dir, recursive = TRUE )
  #######
  ## Let's start to collect the graphs
  ######
  pdf( file=paste0(result_dir,"/" , label4title ,".pdf"),
       width = 10, height = 7)
  meaningful_variables <- setdiff( colnames(melteddf) , c("Unique_ID", "variable", "value") ) # selecting the variables used for the interation
  
  ###########
  # PCAs
  ###########
  for( myfill in meaningful_variables ){
    print(autoplot( prcomp( exp_matrix_T ), data = annotdf, colour=myfill, size=PCA_point_size)+
            ggtitle(paste0(label4title )) )
  }
  
  ##############
  # Boxplots
  ##############
  for( myfill in meaningful_variables ){
    boxp <- ggplot(melteddf, aes_string(y="value", fill = myfill, x = "Unique_ID")) + 
      geom_boxplot( )+ ggtitle(paste( label4title ))  # Density plots
    print(boxp )
  }
  
  ##################
  # Density plots
  ##################
  for( myfill in meaningful_variables ){
    plot1 <- ggplot(melteddf, aes_string(x="value", fill = myfill, y = "Unique_ID")) + 
      geom_density_ridges() + 
      ggtitle(label4title)
    print(plot1) # Density plots
  }
  ####
  ## Make the density plots overlaped like limma does
  ####
  for( myfill in meaningful_variables ){
    plot1 <- ggplot(melteddf ) + 
      geom_density(  aes_string( x="value", group = "Unique_ID" , color=myfill ) ) + 
      ggtitle(label4title)
    print(plot1) # Density plots
  }
  for( myfill in meaningful_variables ){
    colors_from_rainbow <- colors_4_plotDensities(annotdf , myfill )
    plotDensities(exp_matrix, col = colors_from_rainbow )
    #print( plotDensities(exp_matrix, col = colors_from_rainbow ) )  
  }
  #############
  ### heatmap
  #############
  colors_for_heatmap <- colors_4_plotDensities(annotdf ,  group_4_tsnes )
  heatmap( as.matrix(exp_matrix), ColSideColors=colors_for_heatmap , margins=c(30,2), Colv =FALSE)
  heatmap( as.matrix(exp_matrix), ColSideColors=colors_for_heatmap , margins=c(1,2), Colv =FALSE)
  ?heatmap
  ##############
  # Tsne-s
  ##############
  # Tsnes Calculation
  tsne_model_1 = Rtsne(  t(exp_matrix)  , check_duplicates=FALSE, pca=TRUE, perplexity=myperplexity, theta=0.5, dims=2)
  d_tsne_1 = as.data.frame( tsne_model_1$Y )
  
  # Tsne black and white
  d_tsne_1_simplecols <- cbind(d_tsne_1, annotdf[,group_4_tsnes] ) 
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
  
  # Tsne-s colors
  for( myfill in meaningful_variables ){
    d_tsne_1_simplecols <- cbind(d_tsne_1 , annotdf[,myfill] )
    colnames(d_tsne_1_simplecols)[3] <- myfill
    rownames(d_tsne_1_simplecols) <- annotdf$Unique_ID
    mytsneplot_colScanID <- ggplot(d_tsne_1_simplecols, aes_string(x="V1", y="V2" , col=myfill )) +
      geom_point(size=2.5) +
      guides(colour=guide_legend(override.aes=list(size=6))) +
      xlab("") + ylab("") +
      ggtitle( paste("t-SNE",label4title) ) +
      theme_light(base_size=20) +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank()) +
      scale_colour_brewer(palette = "Set1") +
      geom_text(
        label=rownames( d_tsne_1_simplecols ), 
        nudge_x = 0.25, nudge_y = 0.25, 
        check_overlap = T
      )
    print(mytsneplot_colScanID)
  }
  dev.off()  
}
