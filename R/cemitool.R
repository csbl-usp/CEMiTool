#' An S4 class to represent the CEMiTool analysis.
#'
#' @slot expression Gene expression \code{data.frame}.
#' @slot sample_annotation Sample annotation \code{data.frame}.
#' @slot module Genes in modules information \code{data.frame}.
#' @slot enrichment \code{list} with modules enrichment results for sample classes.
#' @slot ora Over-representation analysis results \code{data.frame}.
#' @slot profile_plot list of ggplot graphs with gene expression profile per module.
#' @slot enrichment_plot ggplot graph for enrichment analysis results
#' @slot barplot_ora list of ggplot graphs with over-representation analysis resultsper module
setClass('CEMiTool', slots=list(expression='data.frame',
                                sample_annotation='data.frame',
                                module='data.frame',
                                enrichment='list', # gene set enrichment analysis
                                ora='data.frame',
                                profile_plot='list',
                                enrichment_plot='list',
                                barplot_ora='list',
                                sample_name_column="vector",
                                class_column="vector"))

setMethod("initialize", signature="CEMiTool",
          function(.Object, expression,
                   sample_name_column="SampleName", 
                   class_column="Class", ...){
              .Object@expression <- expression
              .Object@sample_name_column <- sample_name_column
              .Object@class_column <- class_column
              return(.Object)
          })

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
        results <- mod_ora(results, gmt_in=gmt, verbose=verbose)
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
#' @param cem_obj Object of class \code{CEMiTool}
#' 
#' @return number of modules
#'
#' @rdname nmodules
#' @export
setGeneric('nmodules', function(cem_obj) {
    standardGeneric('nmodules')
})

#' @rdname nmodules
setMethod('nmodules', signature(cem_obj='CEMiTool'),
          function(cem_obj) {
              n <- 0
              if(!is.null(cem_obj@module)){
                  n <- length(unique(cem_obj@module$modules))
              } else {
                  message("Run cemitool function to get modules!") 
              }
              return(n)
          }
)

#' Print a cemitool object
#'
#' @export
setMethod('show', signature(object='CEMiTool'),
          function(object) {
              cat("CEMiTool Object\n")
              cat("- Number of modules:", nmodules(object), "\n")
              cat("- Modules: ")
              if(is.null(object@module)){
                  cat("null\n")
              } else {
                  cat("\b\b (data.frame: ", nrow(object@module), "x", ncol(object@module), "): \n", sep="")
                  print(object@module[1:3, ])
              }
              cat("- Gene Set Enrichment Analysis: ")
              if(length(object@enrichment)!=3) {
                  cat("null\n")
              } else {
                  cat("\n    List containing 3 data.frames:\n")
                  cat("        - $ es   : Enrichment Scores ")
                  cat("(", nrow(object@enrichment$es), "x", ncol(object@enrichment$es), ") \n", sep="")
                  cat("        - $ nes  : Normalized Enrichment Scores ")
                  cat("(", nrow(object@enrichment$nes), "x", ncol(object@enrichment$nes), ") \n", sep="")
                  cat("        - $ pval : p-value ")
                  cat("(", nrow(object@enrichment$pval), "x", ncol(object@enrichment$pval), ") \n", sep="")
              }
              cat("- Over Representation Analysis: ")
              if(nrow(object@ora)==0) {
                  cat("null\n")
              } else {
                  cat("\b\b (data.frame: ", nrow(object@ora), "x", ncol(object@ora), "): \n", sep="")
                  print(object@ora[1:3, c('Module', "ID", "p.adjust")])
              }
              cat("- Profile plot: ")
              if(length(object@profile_plot)==0) {
                  cat("null\n")
              } else {
                  cat("ok\n")
              }
              cat("- Enrichment plot: ")
              if(!is.ggplot(object@enrichment_plot)) {
                  cat("null\n")
              } else {
                  cat("ok\n")
              }
              cat("- Barplot of ORA: ")
              if(length(object@barplot_ora)==0) {
                  cat("null\n")
              } else {
                  cat("ok\n")
              }
          }
)


