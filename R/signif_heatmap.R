
#' Draw clustered heatmaps with Significant marker.
#' @param r_data  numeric matrix of the r values to be plotted
#' @param p_data  numeric matrix of the p values to be plotted
#' @param p_level_one  specifies the cutoff of p value to add tag "*"
#' @param p_level_two  specifies the cutoff of p value to add tag "**"
#' @param mcolor  vector of colors used in heatmap
#' @param kmeans_k	 the number of kmeans clusters to make, if we want to aggregate the rows before drawing heatmap. If NA then the rows are not aggregated
#' @param bk a sequence of numbers that covers the range of values in mat and is one element longer than color vector. Used for mapping values to colors. Useful, if needed to map certain values to certain colors, to certain values. If value is NA then the breaks are calculated automatically. When breaks do not cover the range of values, then any value larger than max(breaks) will have the largest color and any value lower than min(breaks) will get the lowest color.
#' @param border_color  color of cell borders on heatmap, use NA if no border should be drawn.
#' @param cellwidth  individual cell width in points. If left as NA, then the values depend on the size of plotting window.
#' @param cellheight 	individual cell height in points. If left as NA, then the values depend on the size of plotting window.
#' @param scale  character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
#' @param cluster_rows 	boolean values determining if rows should be clustered or hclust object,
#' @param cluster_cols 	boolean values determining if columns should be clustered or hclust object.
#' @param clustering_distance_rows  distance measure used in clustering rows. Possible values are "correlation" for Pearson correlation and all the distances supported by dist, such as "euclidean", etc. If the value is none of the above it is assumed that a distance matrix is provided.
#' @param clustering_distance_cols  distance measure used in clustering columns. Possible values the same as for clustering_distance_rows.
#' @param clustering_method 	clustering method used. Accepts the same values as hclust.
#' @param cutree_rows 	number of clusters the rows are divided into, based on the hierarchical clustering (using cutree), if rows are not clustered, the argument is ignored
#' @param cutree_cols  similar to cutree_rows, but for columns
#' @param treeheight_row 	the height of a tree for rows, if these are clustered. Default value 50 points.
#' @param treeheight_col 	the height of a tree for columns, if these are clustered. Default value 50 points.
#' @param legend  logical to determine if legend should be drawn or not.
#' @param legend_breaks  vector of breakpoints for the legend.
#' @param legend_labels  vector of labels for the legend_breaks.
#' @param annotation_row  data frame that specifies the annotations shown on left side of the heatmap. Each row defines the features for a specific row. The rows in the data and in the annotation are matched using corresponding row names. Note that color schemes takes into account if variable is continuous or discrete.
#' @param annotation_col  similar to annotation_row, but for columns.
#' @param annotation  deprecated parameter that currently sets the annotation_col if it is missing
#' @param annotation_colors  list for specifying annotation_row and annotation_col track colors manually. It is possible to define the colors for only some of the features. Check examples for details.
#' @param annotation_legend 	boolean value showing if the legend for annotation tracks should be drawn.
#' @param annotation_names_row  boolean value showing if the names for row annotation tracks should be drawn.
#' @param annotation_names_col 	boolean value showing if the names for column annotation tracks should be drawn.
#' @param drop_levels 	logical to determine if unused levels are also shown in the legend
#' @param show_rownames	 boolean specifying if column names are be shown.
#' @param show_colnames	 boolean specifying if column names are be shown.
#' @param main  the title of the plot
#' @param fontsize  base fontsize for the plot
#' @param fontsize_row  fontsize for rownames (Default: fontsize)
#' @param fontsize_col	fontsize for colnames (Default: fontsize)
#' @param angle_col  angle of the column labels, right now one can choose only from few predefined options (0, 45, 90, 270 and 315)
#' @param number_format  format strings (C printf style) of the numbers shown in cells.
#' @param number_color  color of the text
#' @param fontsize_number  fontsize of the numbers displayed in cells
#' @param gaps_row  vector of row indices that show where to put gaps into heatmap. Used only if the rows are not clustered. See cutree_row to see how to introduce gaps to clustered rows.
#' @param gaps_col  similar to gaps_row, but for columns.
#' @param labels_row  custom labels for rows that are used instead of rownames.
#' @param labels_col  similar to labels_row, but for columns.
#' @param filename  file path where to save the picture. Filetype is decided by the extension in the path. Currently following formats are supported: png, pdf, tiff, bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is calculated so that the plot would fit there, unless specified otherwise.
#' @param width  manual option for determining the output file width in inches.
#' @param height  manual option for determining the output file height in inches.
#' @param silent  do not draw the plot (useful when using the gtable output)
#' @param na_col  specify the color of the NA cell in the matrix.
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer  brewer.pal
#' @export
signif_heatmap<-function(r_data,p_data,p_level_one=0.05,p_level_two=0.01,mcolor = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), kmeans_k = NA, bk = NA, border_color = "grey60",
                         cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE,
                         cluster_cols = TRUE, clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean", clustering_method = "complete",
                         cutree_rows = NA, cutree_cols = NA,
                         treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows,50, 0), treeheight_col = ifelse((class(cluster_cols) == "hclust") ||cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA,
                         legend_labels = NA, annotation_row = NA, annotation_col = NA,
                         annotation = NA, annotation_colors = NA, annotation_legend = TRUE,
                         annotation_names_row = TRUE, annotation_names_col = TRUE,
                         drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA,
                         fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize,
                         angle_col = c("270", "0", "45", "90", "315"),
                         number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8
                         * fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL,
                         labels_col = NULL, filename = NA, width = NA, height = NA,
                         silent = FALSE, na_col = "#DDDDDD"){
  #判断显著性
  if (!is.null(p_data)){
    ssmt <- p_data< p_level_two
    p_data[ssmt] <-'**'
    smt <- p_data >p_level_two& p_data <p_level_one
    p_data[smt] <- '*'
    p_data[!ssmt&!smt]<- ''
  } else {
    p_data <- F
  }
  #可视化
  picture<-pheatmap(r_data,color=mcolor,kmeans_k=kmeans_k,breaks=bk,
                    border_color=border_color,cellwidth=cellwidth,cellheight=cellheight,scale=scale,cluster_rows=cluster_rows,
                    cluster_cols=cluster_cols,clustering_distance_rows=clustering_distance_rows,clustering_distance_cols=clustering_distance_cols,
                    clustering_method=clustering_method,cutree_rows=cutree_rows,cutree_cols=cutree_cols,
                    treeheight_row=treeheight_row,treeheight_col=treeheight_col,legend=legend,legend_breaks=legend_breaks,
                    legend_labels=legend_labels,annotation_row=annotation_row,annotation_col=annotation_col,annotation=annotation,
                    annotation_colors=annotation_colors,annotation_legend=annotation_legend,annotation_names_row=annotation_names_row,
                    annotation_names_col=annotation_names_col,drop_levels=drop_levels,show_rownames=show_rownames,
                    show_colnames=show_colnames,main=main,fontsize=fontsize,fontsize_row=fontsize_row,fontsize_col=fontsize_col,
                    angle_col=angle_col,display_numbers=p_data,number_format=number_format,
                    number_color=number_color,fontsize_number=fontsize_number,gaps_row=gaps_row,gaps_col=gaps_col,
                    labels_row=labels_row,labels_col=labels_col,filename=filename,width=width,height=height,
                    silent=silent,na_col=na_col)
  return(picture)
}
