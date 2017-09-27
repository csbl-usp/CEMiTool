#' @import grid
#' @import gridExtra
#' @import ggplot2
#' @import ggdendro
#' @import gtable
NULL

#' Sample clustering
#' 
#' Creates a dendrogram showing the similarities between samples in the expression data.
#' 
#' @param cem Object of class \code{CEMiTool} or \code{data.frame}.
#' @param col_vector A vector of columns to use for visualizing the clustering. See Details.
#' @param sample_name_column A string specifying the column to be used as sample identification.
#' 		  For CEMiTool objects, this will be the string specified in the sample_name_column slot. 
#' @param filtered Logical. Whether or not to use filtered data for CEMiTool objects (Default: FALSE).
#' @param ... Optional parameters.
#' 
#' @return Object of class \code{CEMiTool} with dendrogram or a plot object. 
#' 
#' @examples 
#' # Get example CEMiTool object
#' data(cem)
#' # Plot sample dendrogram
#' cem <- plot_sample_tree(cem)
#' # Check resulting plot
#' show_plot(cem, "sample_tree")
#'
#' @rdname plot_sample_tree
#' @export
setGeneric('plot_sample_tree', function(cem, ...){
	standardGeneric('plot_sample_tree')
})

#' @rdname plot_sample_tree
setMethod('plot_sample_tree', signature('CEMiTool'),
	function(cem,  col_vector=NULL, sample_name_column=NULL,
			 class_column=NULL, filtered=FALSE){

		expr <- expr_data(cem, filtered=filtered)
	    if(nrow(sample_annotation(cem)) > 0){
		    annot <- sample_annotation(cem)
        	sample_name_column=cem@sample_name_column
		    class_column=cem@class_column
		}else{
			annot <- NULL
		}
    
	    expr_t <- as.data.frame(t(expr))
	    samples <- rownames(expr_t)
	    
	    if(!is.null(annot)){
			rownames(annot) <- annot[, sample_name_column]
		    annot[, sample_name_column] <- NULL
			#annot_class <- annot[, class_column, drop=FALSE]
			annot_rows <- match(samples, rownames(annot))
			if(!is.null(col_vector)){
				if(is.numeric(col_vector) | is.character(col_vector)){
					annot <- annot[, col_vector]
				}else{
				    stop("Accepted classes for col_vector object are 'numeric'.")
				}
			}
			annot <- annot[annot_rows, , drop=FALSE]
								        
			#annot[, "class_code"] <- match(annot$Class, unique(annot[, class_column]))
			#annot_class_code <- annot[, "class_code", drop=FALSE]
			#annot[, "class_code"] <- NULL
								        
		}
		    
		sample_tree <- hclust(dist(expr_t), method = "average")
		#plot(sample_tree)
		colors_samples <- data.frame(class=annot[, class_column],
				                       samples=factor(sample_tree$labels,
		                               levels=sample_tree$labels[sample_tree$order]))
		rownames(colors_samples) <- colors_samples$samples
			    
		lvl <- levels(colors_samples$samples)

		suppressMessages(
			p1 <- ggdendro::ggdendrogram(sample_tree, rotate=FALSE) +
			             scale_x_continuous(expand = c(0, 0.5),
		                 labels = levels(colors_samples$samples),
		                 breaks = 1:length(colors_samples$samples)) +
			             scale_y_continuous(expand = c(0.02, 0))
		)
				    
		p2 <- ggplot2::ggplot(colors_samples, aes(samples, y=1, fill=factor(class))) +
					geom_tile() +
					scale_y_continuous(expand=c(0, 0)) +
					theme(axis.title=element_blank(),
					axis.ticks=element_blank(),
					axis.text=element_blank(),
					legend.position="none")
					    
		annot_num <- Filter(is.numeric, annot)
		order <- match(lvl, rownames(annot_num))
		annot_num <- annot_num[order, ]
							    
		png("NULL")
		gp1 <- ggplot2::ggplotGrob(p1)  
		gp2 <- ggplot2::ggplotGrob(p2)  
								    
		if(ncol(annot_num) > 0){
			df_scaled <- as.data.frame(scale(annot_num))
			df_scaled[, sample_name_column] <- rownames(annot_num)
			df <- reshape2::melt(df_scaled, id.vars=sample_name_column)
			df[, sample_name_column] <- factor(df[, sample_name_column], levels=lvl)

			#custom_pal <- c("#053061", "#2166AC", "#4393C3", "#92C5DE",
			#                "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
			#                "#D6604D", "#B2182B", "#67001F")
			#custom_pal <- colorRampPalette(custom_pal)(200)
															        
			p3 <- ggplot2::ggplot(df, aes_string(x=sample_name_column, y="variable", fill="value")) +
					geom_raster() +
					#scale_y_continuous(expand=c(0, 0)) +
					#scale_fill_gradientn(colors=custom_pal) +
					scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
					#scale_fill_gradientn(colors=blueWhiteRed(100)) +
															        
					theme(axis.title=element_blank(),
					      axis.ticks=element_blank(),
					      axis.text.x=element_blank(),
					      legend.position="none") 
																	        
			gp3 <- ggplot2::ggplotGrob(p3)
			maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5], gp3$widths[2:5])
			gp1$widths[2:5] <- as.list(maxWidth)
			gp2$widths[2:5] <- as.list(maxWidth)
			gp3$widths[2:5] <- as.list(maxWidth)
			invisible(dev.off())
			g <- gridExtra::arrangeGrob(gp1, gp2, gp3, ncol=1,heights=c(2/5, 1/5, 2/5))
			cem@sample_tree_plot <- g
			return(cem)
		}else{
			maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
			gp1$widths[2:5] <- as.list(maxWidth)
			gp2$widths[2:5] <- as.list(maxWidth)
			invisible(dev.off())
			g <- gridExtra::arrangeGrob(gp1, gp2, ncol=1,heights=c(3/5, 2/5))
			cem@sample_tree_plot <- g
			return(cem)
		}
})


