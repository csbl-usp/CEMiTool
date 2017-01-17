# creates a CEMiTool object
setClass('CEMiTool', slots=list(expression='data.frame',
                                sample_annotation='data.frame',
                                module='data.frame',
                                enrichment='list', # gene set enrichment analysis
                                ora='data.frame',
                                profile_plot='list',
                                enrichment_plot='list',
                                barplot_ora='list'))

#' Full gene co-expression analysis
#'
#' Defines co-expression modules and functionally characterizes
#' each one of them.
#'
#' @param exprs Gene expression \code{data.frame}.
#' @param annot Sample annotation \code{data.frame}.
#' @param gmt a list from function prepare_gmt containing the gene sets.
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
#' @return a cemitool object 
#'
#' @examples
#' data(exprs) 
#' cemitool(exprs=exprs)
#'
#' @export
cemitool <- function(exprs, 
                     annot,
                     gmt,
                     cor_method=c('pearson', 'spearman'),
                     merge_similar=TRUE,
                     split_modules=FALSE,
                     ora_pval=0.05,
                     min_ngen=30,
                     diss_thresh=0.8,
                     plot=FALSE,
                     verbose=FALSE)
{
    results <- new('CEMiTool', expression=exprs)
    results@module <- find_modules(results@expression,
                                   cor_method=match.arg(cor_method),
                                   min_ngen=min_ngen,
                                   merge_similar=merge_similar,
                                   diss_thresh=diss_thresh,
                                   verbose=verbose)

    # if user wants splitted modules
    if (split_modules) {
        results <- split_modules(results, min_ngen=min_ngen,
                                 verbose=verbose)
    }

    # if user provides annot file
    if (!missing(annot)) {
        results@sample_annotation <- annot
        #run mod_gsea
        results <- mod_gsea(results, verbose=verbose)
    }

    # if user provides .gmt file
    if (!missing(gmt)) {
        #run mod_ora
        results <- mod_ora(results, gmt=gmt, verbose=verbose)
    }

    # plots all desired charts
    if (plot) {
        results <- plot_profile(results)

        if (!is.null(results@enrichment)) {
            results <- plot_gsea(results)
        }

        if (!is.null(results@ora)) {
            results <- plot_ora(results)
        }
    }
    return(results)
}


#' Get the number of modules on a cemitool object
#'
#' @param x cemitool object
#' 
#' @return number of modules
#'
#' @export
number_modules <- function(x) UseMethod("number_modules")

#' Get the number of modules on a cemitool object
#'
#' @param x cemitool object
#' 
#' @return number of modules
#'
#' @export
number_modules.cemitool <- function(x) {
    n <- 0
    if(!is.null(res@module)){
        n <- length(unique(res@module$modules))
    } else {
        message("Run cemitool function to get modules!") 
    }
    return(n)
}

#' Print a cemitool object
#'
#' @param x cemitool object
#' 
#' @return nothing
#'
#' @export
print.cemitool <- function(x) {
    cat("CEMiTool Object\n")
    cat("- Number of modules:", number_modules(x), "\n")
    cat("- Modules: ")
    if(is.null(x@module)){
        cat("null\n")
    } else {
        cat("\b\b (data.frame: ", nrow(x@module), "x", ncol(x@module), "): \n", sep="")
        print(x@module[1:3, ])
    }
    cat("- Gene Set Enrichment Analysis: ")
    if(is.null(x@enrichment)) {
        cat("null\n")
    } else {
        cat("\n    List containing 3 data.frames:\n")
        cat("        - $ es   : Enrichment Scores ")
        cat("(", nrow(x@enrichment$es), "x", ncol(x@enrichment$es), ") \n", sep="")
        cat("        - $ nes  : Normalized Enrichment Scores ")
        cat("(", nrow(x@enrichment$nes), "x", ncol(x@enrichment$nes), ") \n", sep="")
        cat("        - $ pval : p-value ")
        cat("(", nrow(x@enrichment$pval), "x", ncol(x@enrichment$pval), ") \n", sep="")
    }
    cat("- Over Representation Analysis: ")
    if(is.null(x@ora)) {
        cat("null\n")
    } else {
        cat("\b\b (data.frame: ", nrow(x@ora), "x", ncol(x@ora), "): \n", sep="")
        print(x@ora[1:3, c('Module', "ID", "p.adjust")])
    }
    cat("- Profile plot: ")
    if(is.null(x@profile_plot)) {
        cat("null\n")
    } else {
        cat("ok\n")
    }
    cat("- Enrichment plot: ")
    if(is.null(x@enrichment_plot)) {
        cat("null\n")
    } else {
        cat("ok\n")
    }
    cat("- Barplot of ORA: ")
    if(is.null(x@barplot_ora)) {
        cat("null\n")
    } else {
        cat("ok\n")
    }
}



