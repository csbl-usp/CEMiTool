#' Expression profile visualization
#'
#' Creates a line plot of each gene inside the module through the samples
#'
#' @param exprs Expression data frame or list of expression data frames of each module
#' @param gene_module Two column \code{data.frame}. First column with
#'        gene identifiers and second column with module information
#' @param annot optional. Two column \code{data.frame}. Fist column with
#'        sample identifiers and second column with group/class information
#'
#' @return None
#'
#' @examples
#' plot_profile(exprs)
#'
#' @export
plot_profile <- function(exprs, gene_module, annot, save=FALSE,
                         filename='profile.pdf'){
    modules <- unique(gene_module[, 'modules'])
    plots <- lapply(modules, function(mod){
        # subsets from exprs all genes inside module mod
        genes <- gene_module[gene_module[,'modules']==mod, 'genes']
        mod_exprs <- melt(exprs[genes,])

        g <- ggplot(mod_exprs, aes(x=, y=))

        if (!is.missing(annot)) {
            g <- g + geom_rect()
        }

        g <- g + geom_line(aes(group=genes)) +
                 stat_summary(fun.y=mean, geom='line')

        return(g)
    })
    return(plots)
}



#' Expression profile visualization
#'
#' Creates a line plot of each gene inside the module through the samples
#'
#' @param test
#'
#' @return None
#'
#' @examples
#' plot_graph(test)
#'
#' @export
plot_graph <- function(test){}



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
#' @param test
#'
#' @return None
#'
#' @examples
#' plot_gsea(test)
#'
#' @export
plot_gsea <- function(test){}
