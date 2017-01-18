#' @import ggplot2
NULL

#' Expression profile visualization
#'
#' Creates a line plot of each gene inside the module through the samples
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param ... Optional parameters.
#'
#' @return List with one profile plot per module in gene_module
#'
#' @examples
#' plot_profile(cem_obj)
#'
#' @rdname plot_profile
#' @export
setGeneric('plot_profile', function(cem_obj, ...) {
    standardGeneric('plot_profile')
})

#' @rdname plot_profile
setMethod('plot_profile', signature('CEMiTool'),
          function(cem_obj, order=TRUE) {
              modules <- unique(cem_obj@module[, 'modules'])
              exprs <- cem_obj@expression
              annot <- cem_obj@sample_annotation
              sample_name_column <- cem_obj@sample_name_column
              class_column <- cem_obj@class_column
			  mod_cols <- mod_colors(cem_obj)
              plots <- lapply(modules, function(mod){
                                  # subsets from exprs all genes inside module mod
                                  genes <- cem_obj@module[cem_obj@module[,'modules']==mod, 'genes']
                                  exprs[, 'id'] <- rownames(exprs)
                                  mod_exprs <- melt(exprs[genes,], 'id',
                                                    variable.name='sample',
                                                    value.name='expression')

                                  # initialize plot base layer
                                  g <- ggplot(mod_exprs, aes(x=sample, y=expression))

                                  # adds different background colours if annot is provided
                                  if (nrow(annot)!=0) {
                                      if (order) {
                                          # sorts data.frame by class name
                                          annot <- annot[order(annot[, class_column]),]
                                          annot[, sample_name_column] <- factor(annot[, sample_name_column],
                                                                                levels=annot[, sample_name_column])
                                          mod_exprs[, 'sample'] <- factor(mod_exprs[, 'sample'],
                                                                          levels=annot[, sample_name_column])
                                      }

                                      # y positioning of background tiles
                                      y_pos <- mean(mod_exprs[, 'expression'])

                                      # reinitialize base layer adding background tiles
                                      g <- ggplot(mod_exprs, aes(x=sample, y=expression)) +
                                          geom_tile(data=annot, alpha=0.3, height=Inf,
                                                    aes(x=get(sample_name_column), y=y_pos,
                                                        fill=get(class_column)))
                                  }

                                  # adding lines
                                  g <- g + geom_line(aes(group=id), alpha=0.2, colour=mod_cols[mod]) +
                                      stat_summary(aes(group=1), size=1, fun.y=mean, geom='line')

                                  # custom theme
                                  g <- g + theme(plot.title=element_text(lineheight=0.8,
                                                                         face='bold',
                                                                         colour='black',
                                                                         size=15),
                                                 axis.title=element_text(face='bold',
                                                                         colour='black',
                                                                         size=15),
                                                 axis.text.y=element_text(angle=0,
                                                                          vjust=0.5,
                                                                          size=8),
                                                 axis.text.x=element_text(angle=90,
                                                                          vjust=0.5,
                                                                          size=6),
                                                 panel.grid=element_blank(),
                                                 legend.title=element_blank(),
                                                 legend.text=element_text(size = 8),
                                                 legend.background=element_rect(fill='gray90',
                                                                                size=0.5,
                                                                                linetype='dotted'),
                                                 legend.position='bottom'
                                                 )
                                  # title
                                  g <- g + ggtitle(mod)

                                  return(g)
})
              names(plots) <- modules
              cem_obj@profile_plot <- plots
              return(cem_obj)
          })



#' ORA visualization
#'
#' Creates a bar plot with the results of overenrichment analysis of co-expression modules
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param n number of modules to show
#' @param ... paramaters to plot_ora_single
#'
#' @return a list with ggplot2 object and the number of significant gene sets
#'
#' @examples
#' plot_ora(test)
#'
#' @rdname plot_ora
#' @export
setGeneric('plot_ora', function(cem_obj, ...) {
    standardGeneric('plot_ora')
})

#' @rdname plot_ora
setMethod('plot_ora', signature('CEMiTool'),
          function(cem_obj, n=10, ...){
              ora_splitted <- split(cem_obj@ora, cem_obj@ora$Module)
			  mod_cols <- mod_colors(cem_obj)
              res <- lapply(ora_splitted, function(x){
                                plot_ora_single(head(x, n=n),
                                                graph_color=mod_cols[unique(x$Module)],
                                                title=unique(x$Module),
                                                ...)
              })
              cem_obj@barplot_ora <- res
              return(cem_obj)
          }
)


#' ORA visualization for one module
#'
#' @param es a data.frame from ora function containing only one module
#' @param ordr_by column to order the data.frame
#' @param max_length max length of a gene set name
#' @param pv_cut p-value cuttoff
#' @param graph_color color of bars
#' @param title title of the graph
#'
#' @return a list with ggplot2 object and the number of significant gene sets
#'
#' @examples
#' plot_ora_single(test)
#'
plot_ora_single <- function(es, ordr_by='p.adjust', max_length=50, pv_cut=0.01,
                            graph_color="#4169E1", title="Over Representation Analysis"){
    
    es[, "GeneSet"] <- es[, "ID"]
    
    # limits name length
    ovf_rows <- which(nchar(es[, "GeneSet"]) > max_length) # overflow
    es[ovf_rows, "GeneSet"] <-  paste0(strtrim(es[ovf_rows, "GeneSet"], max_length), "...")
    es[, "GeneSet"] <- stringr::str_wrap(es[, "GeneSet"], width = 20)
    
    # order bars
    lvls <- es[order(es[, ordr_by], decreasing=TRUE), "GeneSet"]
    es[, "GeneSet"] <- factor(es[, "GeneSet"], levels=lvls)
    
    es[, "alpha"] <- 1
    es[es[, ordr_by] > pv_cut, "alpha"] <- 0
    
    # Avoid 0's
    es[es[, ordr_by] > 0.8, ordr_by] <- 0.8
    my_squish <- function(...){
        return(scales::squish(..., only.finite=FALSE))
    }
    
    # plot
    y_axis <- paste('-log10(', ordr_by, ')')
    pl <- ggplot(es, aes_string(x="GeneSet", y=y_axis, alpha="alpha", fill=y_axis)) + 
        geom_bar(stat="identity") +
        theme(axis.text=element_text(size=8), legend.title=element_blank()) +
        coord_flip() + 
        scale_alpha(range=c(0.4, 1), guide="none") +
        labs(y="-log10(adjusted p-value)", title=title, x="") +
        geom_hline(yintercept=-log10(pv_cut), colour="grey", linetype="longdash") + 
        scale_fill_gradient(low="gray", high=graph_color, limits=c(2, 5), oob=my_squish)
    res <- list('pl'=pl, numsig=sum(es[, ordr_by] < pv_cut, na.rm=TRUE))
    return(res)
}



#' GSEA visualization
#'
#' Creates a heatmap with the results of gene set enrichment analysis (GSEA) of co-expression modules
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param pv_cut P-value cut-off. Default \code{0.05}
#'
#' @return None
#'
#' @examples
#' plot_gsea(enrichment)
#'
#' @rdname plot_gsea
#' @export
setGeneric('plot_gsea', function(cem_obj, ...) {
    standardGeneric('plot_gsea')
})

#' @rdname plot_gsea
setMethod('plot_gsea', signature('CEMiTool'),
          function(cem_obj, pv_cut=0.05) {
              stats <- names(cem_obj@enrichment)
              enrichment <- lapply(cem_obj@enrichment, function(stat){
                                       stat[is.na(stat)] <- 0
                                       rownames(stat) <- stat[,1]
                                       stat[,1] <- NULL
                                       return(stat)
              })
              names(enrichment) <- stats

              pval <- enrichment[['pval']]
              nes <- enrichment[['nes']]

              pval <- pval[rowSums(pval < pv_cut) >= 1,]
              nes <- nes[rownames(pval),]

              # check if there is any signif. module
              if(nrow(nes) < 0){
                  stop("No significant modules found!")
              }

              custom_pal <- c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                              "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                              "#D6604D", "#B2182B", "#67001F")
              custom_pal <- colorRampPalette(custom_pal)(200)

              nes <- as.matrix(nes)
              pval <- as.matrix(pval)
              nes[which(pval > pv_cut, arr.ind=T)] <- 0

              row_order <- rownames(nes)[hclust(dist(nes))$order]

              nes_melted <- reshape2::melt(nes)
              colnames(nes_melted) <- c("Module", "Class", "NES")
              nes_melted$Module <- factor(nes_melted$Module, levels=row_order)
              max_abs_nes <- max(abs(nes_melted$NES))
              res <- ggplot(nes_melted, aes(x=Class, y=Module, size=abs(NES), fill=NES)) + 
                  geom_point(color = "white", shape=21) +
                      scale_fill_gradientn(colours=custom_pal, space = "Lab", 
                                           limits=c(-max_abs_nes, max_abs_nes)) +
          scale_size(range=c(0,40)) +
          guides(size="none") +
          theme_minimal() +
          theme(panel.grid.major = element_blank()) +
          scale_x_discrete(position = "top")
      cem_obj@enrichment_plot <- as.list(res)

      return(cem_obj)
          }
          )

