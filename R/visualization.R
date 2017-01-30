#' @import ggplot2
#' @import ggnetwork
NULL

#' Expression profile visualization
#'
#' Creates a line plot of each gene inside the module through the samples
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param order Logical. If TRUE, sorts samples by class.
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with profile plots 
#'
#' @examples
#' plot_profile(cem)
#'
#' @rdname plot_profile
#' @export
setGeneric('plot_profile', function(cem, ...) {
    standardGeneric('plot_profile')
})

#' @rdname plot_profile
setMethod('plot_profile', signature('CEMiTool'),
          function(cem, order=TRUE) {
              modules <- unique(cem@module[, 'modules'])
              modules <- modules[order(as.numeric(stringr::str_extract(modules, "\\d+")))]
              expr <- expr_data(cem)
              annot <- cem@sample_annotation
              sample_name_column <- cem@sample_name_column
              class_column <- cem@class_column
              mod_cols <- mod_colors(cem)
              plots <- lapply(modules, function(mod){
                                  # subsets from expr all genes inside module mod
                                  genes <- cem@module[cem@module[,'modules']==mod, 'genes']
                                  expr[, 'id'] <- rownames(expr)
                                  mod_expr <- melt(expr[genes,], 'id',
                                                    variable.name='sample',
                                                    value.name='expression')

                                  # initialize plot base layer
                                  g <- ggplot(mod_expr, aes(x=sample, y=expression))

                                  # adds different background colours if annot is provided
                                  if (nrow(annot)!=0) {
                                      if (order) {
                                          # sorts data.frame by class name
                                          annot <- annot[order(annot[, class_column]),]
                                          annot[, sample_name_column] <- factor(annot[, sample_name_column],
                                                                                levels=annot[, sample_name_column])
                                          mod_expr[, 'sample'] <- factor(mod_expr[, 'sample'],
                                                                          levels=annot[, sample_name_column])
                                      }

                                      # y positioning of background tiles
                                      y_pos <- mean(mod_expr[, 'expression'])

                                      # reinitialize base layer adding background tiles
                                      g <- ggplot(mod_expr, aes(x=sample, y=expression)) +
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
              cem@profile_plot <- plots
              return(cem)
          })



#' ORA visualization
#'
#' Creates a bar plot with the results of overenrichment analysis of co-expression modules
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param n number of modules to show
#' @param ... paramaters to plot_ora_single
#'
#' @return Object of class \code{CEMiTool} with ORA plots
#'
#' @examples
#' plot_ora(test)
#'
#' @rdname plot_ora
#' @export
setGeneric('plot_ora', function(cem, ...) {
    standardGeneric('plot_ora')
})

#' @rdname plot_ora
setMethod('plot_ora', signature('CEMiTool'),
          function(cem, n=10, ...){
              ora_splitted <- split(cem@ora, cem@ora$Module)
              mod_cols <- mod_colors(cem)
              res <- lapply(ora_splitted, function(x){
                                plot_ora_single(head(x, n=n),
                                                graph_color=mod_cols[unique(x$Module)],
                                                title=unique(x$Module),
                                                ...)
              })
              cem@barplot_ora <- res
              return(cem)
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
#' @param cem Object of class \code{CEMiTool}.
#' @param pv_cut P-value cut-off. Default \code{0.05}
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with GSEA plots
#'
#' @examples
#' plot_gsea(enrichment)
#'
#' @rdname plot_gsea
#' @export
setGeneric('plot_gsea', function(cem, ...) {
    standardGeneric('plot_gsea')
})

#' @rdname plot_gsea
setMethod('plot_gsea', signature('CEMiTool'),
          function(cem, pv_cut=0.05) {
              stats <- names(cem@enrichment)
              enrichment <- lapply(cem@enrichment, function(stat){
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
                  scale_size(range=c(0,30)) +
                  guides(size="none") +
                  theme_minimal() +
                  theme(panel.grid.major = element_blank()) +
                  scale_x_discrete(position = "top")
              cem@enrichment_plot <- res

              return(cem)
          }
          )


#' Network visualization
#'
#' Creates a graph based on interactions provided
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param n number of nodes to label
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with profile plots 
#'
#' @examples
#' plot_interactions(cem)
#'
#' @rdname plot_interactions
#' @export
setGeneric('plot_interactions', function(cem, ...) {
    standardGeneric('plot_interactions')
})

#' @rdname plot_interactions
setMethod('plot_interactions', signature('CEMiTool'),
          function(cem, n=10, ...) {
              mod_cols <- mod_colors(cem)
              mod_names <- names(cem@interactions)
              hubs <- get_hubs(cem)
              res <- lapply(mod_names, function(name){
                                plot_interaction(ig_obj=cem@interactions[[name]], 
                                                 n=n, color=mod_cols[name], title=name,
                                                 coexp_hubs=hubs[[name]])
                                       })
              names(res) <- mod_names
              cem@interaction_plot <- res
              return(cem)
          })

plot_interaction <- function(ig_obj, n, color, title, coexp_hubs){
    degrees <- igraph::degree(ig_obj, normalized=F)
    max_n <- min(n, length(degrees))
    ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
    net <- ggnetwork(ig_obj)
    net[, "shouldLabel"] <- FALSE
    net[, "Hub"] <- ""
    int_hubs <- names(sort(degrees, decreasing=T))[1:max_n]
    int_bool <- net[, "vertex.names"] %in% int_hubs
    net[which(int_bool), "Hub"] <- "Interaction"
    sel_vertex <- int_hubs
    if(!missing(coexp_hubs)){
        coexp_bool <- net[, "vertex.names"] %in% coexp_hubs
        coexp_and_int <- coexp_bool & int_bool
        net[which(coexp_bool), "Hub"] <- "Co-expression"
        net[which(coexp_and_int), "Hub"] <- "Co-expression + Interaction"
        sel_vertex <- c(sel_vertex, coexp_hubs)
    }
    net[which(net[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <- TRUE
    pl <- ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(color = "#DDDDDD", alpha=0.5) +
        geom_nodes(aes(alpha=degree, size=degree), color=color) +
        geom_nodelabel_repel(aes(label = vertex.names, color=Hub),
                             box.padding = unit(1, "lines"),
                             data = function(x) { x[ x$shouldLabel, ]}) + 
        scale_colour_manual(values=c("Co-expression" = "#005E87",
                                     "Interaction" = "#540814",
                                     "Co-expression + Interaction" = "#736E0B")) +
        labs(title=title) +
        theme_blank()
    return(pl)
}
