#' @import ggplot2
NULL
#' Expression profile visualization
#'
#' Creates a line plot of each gene inside the module through the samples
#'
#' @param exprs Expression data frame or list of expression data frames of each 
#'        module.
#' @param gene_module Two column \code{data.frame}. First column with
#'        gene identifiers and second column with module information
#' @param annot Optional. Two column \code{data.frame}. Fist column with
#'        sample identifiers and second column with group/class information
#' @param sample_col character string with the name of samples column in the
#'        annotation \code{data.frame}. If annot is given, this parameters
#'        required.
#' @param class_col character string with the name of classes column in the
#'        annotation \code{data.frame}. If annot is given, this parameters
#'        required.
#'
#' @return List with one profile plot per module in gene_module
#'
#' @examples
#' plot_profile(exprs, gene_module)
#'
#' @export
plot_profile <- function(exprs, gene_module, annot=NULL, sample_col=NULL,
                         class_col=NULL, order=TRUE, filename='profile.pdf')
{
    modules <- unique(gene_module[, 'modules'])
    plots <- lapply(modules, function(mod){
        # subsets from exprs all genes inside module mod
        genes <- gene_module[gene_module[,'modules']==mod, 'genes']
        exprs[, 'id'] <- rownames(exprs)
        mod_exprs <- melt(exprs[genes,], 'id',
                          variable.name='sample',
                          value.name='expression')

        # initialize plot base layer
        g <- ggplot2::ggplot(mod_exprs, aes(x=sample, y=expression))

        # adds different background colours if annot is provided
        if (!is.null(annot)) {
            if (is.null(sample_col)) {
                stop('Must provide sample column name')
            }
            if (is.null(class_col)) {
                stop('Must provide class column name')
            }

            if (order) {
                # sorts data.frame by class name
                annot <- annot[order(annot[, class_col]),]
                annot[, sample_col] <- factor(annot[, sample_col],
                                            levels=annot[, sample_col])
                mod_exprs[, 'sample'] <- factor(mod_exprs[, 'sample'],
                                                levels=annot[, sample_col])
            }

            # y positioning of background tiles
            y_pos <- mean(mod_exprs[, 'expression'])

            # reinitialize base layer adding background tiles
            g <- ggplot2::ggplot(mod_exprs, aes(x=sample, y=expression)) +
                        geom_tile(data=annot, alpha=0.3, height=Inf,
                                  aes(x=get(sample_col), y=y_pos,
                                      fill=get(class_col)))
        }

        # adding lines
        g <- g + ggplot2::geom_line(aes(group=id), alpha=0.2) +
                 ggplot2::stat_summary(aes(group=1),size=1,fun.y=mean,geom='line')

        # custom theme
        g <- g + ggplot2::theme(plot.title=ggplot2::element_text(lineheight=0.8,
                                                                 face='bold',
                                                                 colour='black',
                                                                 size=15),
                                axis.title=ggplot2::element_text(face='bold',
                                                                 colour='black',
                                                                 size=15),
                                axis.text.y=ggplot2::element_text(angle=0,
                                                                  vjust=0.5,
                                                                  size=8),
                                axis.text.x=ggplot2::element_text(angle=90,
                                                                  vjust=0.5,
                                                                  size=6),
                       panel.grid=ggplot2::element_blank(),
                       legend.title=ggplot2::element_blank(),
                       legend.text=ggplot2::element_text(size = 8),
                       legend.background=ggplot2::element_rect(fill='gray90',
                                                               size=0.5,
                                                               linetype='dotted'),
                       legend.position='bottom'
                       )
        # title
        g <- g + ggplot2::ggtitle(mod)

        return(g)
    })
    names(plots) <- modules
    return(plots)
}



#' ORA visualization
#'
#' Creates a bar plot with the results of overenrichment analysis of co-expression modules
#'
#' @param ora_res a data.frame from ora function 
#' @param n number of modules to show
#' @param ... paramaters to plot_ora_single
#'
#' @return a list with ggplot2 object and the number of significant gene sets
#'
#' @examples
#' plot_ora(test)
#'
#' @export
plot_ora <- function(ora_res, n=10, ...){
    ora_splitted <- split(ora_res, ora_res$Module)
    res <- lapply(ora_splitted, function(x){
        plot_ora_single(head(x, n=n), 
                        graph_color=sub("\\.\\S$", "", unique(x$Module)),
                        title=unique(x$Module),
                        ...)
    })
    return(res)
}


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
#' @export
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
    pl <- ggplot2::ggplot(es, ggplot2::aes_string(x="GeneSet", y=y_axis, alpha="alpha", fill=y_axis)) + 
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme(axis.text=ggplot2::element_text(size=8), legend.title=ggplot2::element_blank()) +
        ggplot2::coord_flip() + 
        ggplot2::scale_alpha(range=c(0.4, 1), guide="none") +
        ggplot2::labs(y="-log10(adjusted p-value)", title=title, x="") +
        ggplot2::geom_hline(yintercept=-log10(pv_cut), colour="grey", linetype="longdash") + 
        ggplot2::scale_fill_gradient(low="gray", high=graph_color, limits=c(2, 5), oob=my_squish)
    res <- list('pl'=pl, numsig=sum(es[, ordr_by] < pv_cut, na.rm=TRUE))
    return(res)
}



#' GSEA visualization
#'
#' Creates a heatmap with the results of gene set enrichment analysis (GSEA) of co-expression modules
#'
#' @param enrichment Output of mod_gsea function.
#' @param pv_cut P-value cut-off. Default \code{0.05}
#'
#' @return None
#'
#' @examples
#' plot_gsea(enrichment)
#'
#' @export
plot_gsea <- function(enrichment, pv_cut=0.05)
{
    stats <- names(enrichment)
    enrichment <- lapply(enrichment, function(stat){
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
    res <- ggplot2::ggplot(nes_melted, ggplot2::aes(x=Class, y=Module, size=abs(NES), fill=NES)) + 
        ggplot2::geom_point(color = "white", shape=21) +
        ggplot2::scale_fill_gradientn(colours=custom_pal, space = "Lab", 
                                      limits=c(-max_abs_nes, max_abs_nes)) +
        ggplot2::scale_size(range=c(0,40)) +
        ggplot2::guides(size="none") +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank()) +
        ggplot2::scale_x_discrete(position = "top")

    return(res)
}

