#' @import grid
#' @import gridExtra
#' @import ggplot2
#' @import ggdendro
#' @import gtable
#' @import ggpmisc
#' @import ggthemes
#' @importFrom rmarkdown render
#' @import knitr
#' @importFrom DT datatable
#' @import htmltools
NULL

#' Sample clustering
#' 
#' Creates a dendrogram showing the similarities between samples in the expression data.
#' 
#' @param cem Object of class \code{CEMiTool} or \code{data.frame}.
#' @param col_vector A vector of columns to use for visualizing the clustering. See Details.
#' @param sample_name_column A string specifying the column to be used as sample identification.
#' 		  For CEMiTool objects, this will be the string specified in the sample_name_column slot.
#' @param class_column A string specifying the column to be used as sample group identification.
#' 		  For CEMiTool objects, this will be the string specified in the class_column slot.
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
	function(cem, col_vector=NULL, sample_name_column=NULL,
			 class_column=NULL, filtered=FALSE){

		expr <- expr_data(cem, filtered=filtered)
		if(nrow(expr) == 0){
		    stop("CEMiTool object has no expression file!")
	    }
		vars <- mget(ls())
		vars$expr <- NULL
		cem <- get_args(cem, vars)
	    if(nrow(sample_annotation(cem)) > 0){
		    annot <- sample_annotation(cem)
        	sample_name_column <- cem@sample_name_column
		    class_column <- cem@class_column
		}else{
			annot <- NULL
		}
    
	    expr_t <- as.data.frame(t(expr))
	    samples <- rownames(expr_t)
		sample_tree <- hclust(dist(expr_t), method = "average")

	    if(!is.null(annot)){
			rownames(annot) <- annot[, sample_name_column]
		    annot[, sample_name_column] <- NULL
			annot_rows <- match(samples, rownames(annot))
			if(!is.null(col_vector)){
				if(is.numeric(col_vector) | is.character(col_vector)){
					annot <- annot[, col_vector]
				}else{
				    stop("Accepted classes for col_vector object are 'numeric'.")
				}
			}
			annot <- annot[annot_rows, , drop=FALSE]
								        
			colors_samples <- data.frame(class=annot[, class_column],
				                       samples=factor(sample_tree$labels,
		                               			 levels=sample_tree$labels[sample_tree$order]))
		}else{
			colors_samples <- data.frame(samples=factor(sample_tree$labels,
		                               			   levels=sample_tree$labels[sample_tree$order]))
		}
		rownames(colors_samples) <- colors_samples$samples
			    
		lvl <- levels(colors_samples$samples)

		suppressMessages(
			p1 <- ggdendro::ggdendrogram(sample_tree, rotate=FALSE) +
			             scale_x_continuous(expand = c(0, 0.5),
		                 labels = levels(colors_samples$samples),
		                 breaks = 1:length(colors_samples$samples)) +
			             scale_y_continuous(expand = c(0.02, 0))
		)

		png("NULL")
        gp1 <- ggplot2::ggplotGrob(p1)

		if(!is.null(annot)){		    
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
								    
			#png("NULL")
			#gp1 <- ggplot2::ggplotGrob(p1)  
			gp2 <- ggplot2::ggplotGrob(p2)  
									    
			if(ncol(annot_num) > 0){
				df_scaled <- as.data.frame(scale(annot_num))
				df_scaled[, sample_name_column] <- rownames(annot_num)
				df <- reshape2::melt(df_scaled, id.vars=sample_name_column)
				df[, sample_name_column] <- factor(df[, sample_name_column], levels=lvl)
	
																        
				p3 <- ggplot2::ggplot(df, aes_string(x=sample_name_column, y="variable", fill="value")) +
						geom_raster() +
						scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
																        
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
		}else{
			invisible(dev.off())
			g <- gridExtra::arrangeGrob(gp1, ncol=1)
			cem@sample_tree_plot <- g
			return(cem)
		}
})

#' Plot mean and variance
#'
#' This plot returns a scatterplot of the mean by the variance
#' of gene expression. A linear relationship between these values for
#' RNAseq data suggest that an appropriate transformation such as the
#' Variance Stabilizing Transformation should be applied.
#' 
#' @param cem Object of class \code{CEMiTool}
#' @param filtered Logical. Whether or not to use filtered data for CEMiTool objects (Default: FALSE).
#' @param ... Optional parameters
#'
#' @return Object of class \code{CEMiTool} containing a mean and variance plot
#' 
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Plot mean and variance plot
#' cem <- plot_mean_var(cem)
#' # Check results
#' show_plot(cem, 'mean_var')
#'
#' @rdname plot_mean_var
#' @export
setGeneric('plot_mean_var', function(cem, ...){
	standardGeneric('plot_mean_var')
})

#' @rdname plot_mean_var
setMethod('plot_mean_var', signature('CEMiTool'),
	function(cem, filtered=FALSE){
		expr <- expr_data(cem, filtered=filtered)
		if(nrow(expr) == 0){
		    stop("CEMiTool object has no expression file!")
	    }
		vars <- mget(ls())
        vars$expr <- NULL
        cem <- get_args(cem, vars=vars)
	    
	    expr_mean <- apply(expr, 1, mean)
		expr_var <- apply(expr, 1, var)
	
		mean_var <- data.frame(Mean=expr_mean, Variance=expr_var)
		log_mean_var <- as.data.frame(apply(mean_var, 2, log10))
		my_formula <- y ~ x
				    
		pl <- ggplot(log_mean_var, aes(x=Mean, y=Variance)) + 
				geom_point() + 
				geom_smooth(method="lm", se=FALSE, color="red", formula=my_formula)+
				ggpmisc::stat_poly_eq(formula=my_formula, 
					aes(label=paste(..eq.label.., ..rr.label.., sep="*plain(\",\")~")), 
					parse=TRUE) +
	        labs(x = "Mean Expression (log10)", y="Variance (log10)") +
			ggthemes::theme_gdocs() +
			theme(rect=element_blank(),
				  axis.title.x = element_text(face="bold", size=12),
			      axis.title.y = element_text(face="bold", size=12))

		cem@mean_var_plot <- pl
		return(cem)
})

#' Plot histogram
#'
#' This function plots a histogram of the distribution of gene expression, to 
#' help assess the normality of the data. 
#' 
#' @param cem Object of class \code{CEMiTool}
#' @param filtered Logical. Whether or not to use filtered data for CEMiTool objects (Default: FALSE).
#' @param ... Optional parameters
#' 
#' @return Object of class \code{CEMiTool} containing expression histogram
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Plot histogram
#' cem <- plot_hist(cem)
#' # Check results
#' show_plot(cem, "hist")
#'
#' @rdname plot_hist
#' @export

setGeneric('plot_hist', function(cem, ...){
	standardGeneric('plot_hist')
})

#' @rdname plot_hist
setMethod('plot_hist', signature('CEMiTool'),
	function(cem, filtered=FALSE){
		expr <- expr_data(cem, filtered=filtered)
		if(nrow(expr) == 0){
			stop("CEMiTool object has no expression file!")
	    }
        vars <- mget(ls())
        vars$expr <- NULL
        cem <- get_args(cem, vars=vars)
		measures <- as.data.frame(as.vector(as.matrix(expr)))
		names(measures) <- "data"
		minExp <- round(min(measures, na.rm=TRUE)-0.5,digits=0)
		maxExp <- round(max(measures, na.rm=TRUE)+0.5, digits=0)
		delta <- (maxExp -minExp)/100

		pl <- ggplot(measures, aes(measures$data)) + 
				geom_histogram(breaks=seq(minExp,maxExp , by=delta), 
							   col="lightgrey",
							   fill="#4A7CB2") +
		     	labs(x="Measures", y="Count") +
			    ggthemes::theme_gdocs() +
				theme(rect=element_blank(),
					  axis.title.x = element_text(face="bold", size=12),
					  axis.title.y = element_text(face="bold", size=12), 
					  panel.grid=element_blank())

		cem@hist_plot <- pl
		return(cem)
})

#' Plot quantile-quantile plot
#'
#' This function creates a normal QQ plot of the expression values.
#'
#' @param cem Object of class \code{CEMiTool}
#' @param filtered Logical. Whether or not to use filtered data for CEMiTool objects (Default: FALSE).
#' @param ... Optional parameters
#' 
#' @return Object of class \code{CEMiTool} containing qqplot
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Plot quantile-quantile plot
#' cem <- plot_qq(cem)
#' # Check results
#' show_plot(cem, 'qq')
#' 
#' @rdname plot_qq
#' @export

setGeneric('plot_qq', function(cem, ...){
	standardGeneric('plot_qq')
})

#' @rdname plot_qq
setMethod('plot_qq', signature('CEMiTool'),
	function(cem, filtered=FALSE){
		expr <- expr_data(cem, filtered=filtered)
	        if(nrow(expr) == 0){
			    stop("CEMiTool object has no expression file!")
	    }
        vars <- mget(ls())
        vars$expr <- NULL
        cem <- get_args(cem, vars=vars)
		measures <- as.data.frame(as.vector(as.matrix(expr)))
		names(measures) <- "data"

		pl <- ggplot(measures, aes(sample = data)) + 
		    stat_qq() + 
		    stat_qq_line() +
		    ggthemes::theme_gdocs() +
		    theme(rect=element_blank(),
				  axis.title.x = element_text(face="bold", size=12),
	              axis.title.y = element_text(face="bold", size=12), 
		          panel.grid=element_blank())

		cem@qq_plot <- pl
		return(cem)
})

#' Diagnostic report
#'
#' Creates report for identifying potential problems with data.
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param directory Directory name for results.
#' @param title Character string with the title of the report.
#' @param force If the directory exists, execution will not stop. 
#' @param ... parameters to rmarkdown::render
#' 
#' @return An HTML file with an interactive diagnostic report.
#'
#' @rdname diagnostic_report
#' @export
setGeneric('diagnostic_report', function(cem, ...) {
	standardGeneric('diagnostic_report')
})

#' @rdname diagnostic_report
setMethod('diagnostic_report', signature('CEMiTool'),
	function(cem, title="Diagnostics", directory="./Reports/Diagnostics", force=FALSE, ...) {
		if(dir.exists(directory)){
			if(!force){
				stop("Stopping analysis: ", directory, " already exists! Use force=TRUE to overwrite.")
			}
		}else{
			dir.create(directory, recursive=TRUE)
		}
		rmd <- system.file("diagnostics", "diagnostics.Rmd", package = "CEMiTool")
		rmarkdown::render(rmd, output_dir=directory, intermediates_dir=directory, ...)
})







