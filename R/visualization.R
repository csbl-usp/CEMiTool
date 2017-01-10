suppressPackageStartupMessages({
    require(igraph)
    require(ggplot2)
    require(ggnetwork)
    require(corrplot)
})

#' Expression profile visualization
#'
#' Creates a line plot of each gene inside the module through the samples
#'
#' @param exprs Expression data frame or list of expression data frames of each module
#' @param gene_module Two column \code{data.frame}. First column with
#'        gene identifiers and second column with module information
#' @param annot optional. Two column \code{data.frame}. Fist column with
#'        sample identifiers and second column with group/class information
#' @param sample_col character string with the name of samples column in the
#'        annotation data.frame
#' @param class_col character string with the name of classes column in the
#'        annotation data.frame
#'
#' @return List with one profile plot per module in gene_module
#'
#' @examples
#' plot_profile(exprs, coex)
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
        g <- ggplot(mod_exprs, aes(x=sample, y=expression))

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
            g <- ggplot(mod_exprs, aes(x=sample, y=expression)) +
                        geom_tile(data=annot, alpha=0.3, height=Inf,
                                  aes(x=get(sample_col), y=y_pos,
                                      fill=get(class_col)))
        }

        # adding lines
        g <- g + geom_line(aes(group=id), alpha=0.2) +
                 stat_summary(aes(group=1),size=1, fun.y=mean, geom='line')

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
    return(plots)
}



#' Expression profile visualization
#'
#' Creates a line plot of each gene inside the module through the samples
#'
#' @param edgelist Two column \code{data.frame} with genes interaction data
#'
#' @return List with one network plot per module in gene_module
#'
#' @examples
#' plot_graph(edgelist)
#'
#' @export
plot_graph <- function(edgelist){}



#' ORA visualization
#'
#' Creates a bar plot with the results of overenrichment analysis of co-expression modules
#'
#' @param test
#'
#' @return None
#'
#' @examples
#' plot_ora(test)
#'
#' @export
plot_ora <- function(test){}



#' GSEA visualization
#'
#' Creates a heatmap with the results of gene set enrichment analysis (GSEA) of co-expression modules
#'
#' @param enrichment
#' @param pv_cut
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

    order <- rownames(nes)[hclust(dist(nes))$order]

    corrplot(nes[order, ], p.mat=pval[order, ], col=custom_pal,
             is.corr=FALSE, addgrid.col="white", insig="blank",
             pch.cex=0.5, pch.col="black", tl.col="black", tl.cex=0.5,
             cl.cex=0.4, cl.ratio=0.5, cl.pos="r", cl.align.text="l",
             mar=c(0,0,0,0), sig.level=pv_cut)
}
