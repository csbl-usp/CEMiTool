#' Full gene co-expression analysis
#'
#' Defines co-expression modules and functionally characterizes
#' each one of them.
#'
#' @param exprs Gene expression \code{data.frame}.
#' @param annot Sample annotation \code{data.frame}.
#' @param gmt A character string with name of the Gene set file in GMT format.
#' @param cor_method A character string indicating which correlation coefficient is
#'        to be computed. One of \code{"pearson"} or \code{"spearman"}.
#'        Default \code{"pearson"}.
#' @param merge_similar Logical. If \code{TRUE}, merge similar modules.
#' @param split_modules Logical. If \code{TRUE}, splits modules by correlation sign.
#' @param ora_pval P-value for overrepresentation analysis. Default \code{0.05}.
#' @param min_ngen Minimal number of genes per submodule. Default \code{30}.
#' @param diss_thresh Module merging correlation threshold for eigengene similarity.
#'        Default \code{0.8}.
#' @param plot Logical. If \code{TRUE}, plots all figures.
#' @param verbose Logical. If \code{TRUE}, reports analysis steps.
#'
#' @return just god knows 
#'
#' @examples
#' cemitool(exprs=expression.df)
#'
#' @export
cemitool <- function(exprs, 
                     annot,
                     gmt,
                     cor_method=c('pearson', 'spearman'),
                     merge_similar=TRUE,
                     split_mods=FALSE,
                     ora_pval=0.05,
                     min_ngen=30,
                     diss_thresh=0.8,
                     plot=FALSE,
                     verbose=FALSE)
{
    gene_module <- find_modules(exprs,
                                cor_method=match.arg(cor_method),
                                min_ngen=min_ngen,
                                merge_similar=merge_similar,
                                diss_thresh=diss_thresh,
                                verbose=verbose)
    
    # if user wants splitted modules
    if (split_mods) {
        gene_module <- split_modules(exprs=exprs, gene_module=gene_module,
                                     min_ngen=min_ngen,
                                     verbose=verbose)
    }

    # if user provides annot file
    if (!is.null(annot)) {
        #run mod_gsea
        gsea <- mod_gsea(exprs=exprs, gene_module=gene_module,
                         annot=annot, verbose=verbose)
    }

    # if user provides .gmt file
    if (!is.null(gmt)) {
        #run mod_ora
        ora <- mod_ora(gene_module=gene_module, gmt=gmt, verbose=verbose)
    }

    # plots all desired charts
    if (plot) {
        profiles <- plot_profile(exprs, gene_module)
        if (exists('gsea')) {
            pdf('modules_enrichment.pdf')
            print(plot_gsea(gsea))
        }

        if (exists('ora')) {
            
        }


    }

}


