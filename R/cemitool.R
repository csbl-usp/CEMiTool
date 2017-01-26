#' @importFrom grDevices rainbow
#' @importFrom utils head
#' @importFrom methods new 'slot<-' show

setOldClass('gg')
setOldClass('ggplot')
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
                                selected_genes='vector',
                                module='data.frame',
                                enrichment='list', # gene set enrichment analysis
                                ora='data.frame',
                                interactions='list',
                                interaction_plot='list',
                                profile_plot='list',
                                enrichment_plot='gg',
                                barplot_ora='list',
                                sample_name_column='vector',
                                class_column='vector',
                                mod_colors='character',
                                parameters='list',
                                adjacency='matrix'))

setMethod("initialize", signature="CEMiTool",
          function(.Object, expression,
                   sample_name_column="SampleName", 
                   class_column="Class", ...){
              .Object@sample_name_column <- sample_name_column
              .Object@class_column <- class_column
              arguments <- list(...)
              for( arg_name in names(arguments) ) {
                  slot(.Object, arg_name) <- arguments[[arg_name]]
              }
              if (!missing(expression)) {
                  slot(.Object, 'expression') <- expression
                  slot(.Object, 'selected_genes') <- rownames(expression) 
              }
              return(.Object)
          })


#' Retrieve and set expression attribute
#'
#' @param cem_obj Object of class \code{CEMiTool}
#' @param value Object of class \code{data.frame} with gene
#'        expression data
#' @param filtered logical. If TRUE retrieves filtered expression data
#' 
#' @return Object of class \code{data.frame} with gene expression data
#'
#' @rdname expr_data 
#' @export
setGeneric("expr_data", function(cem_obj, ...) {
            standardGeneric("expr_data")
          })

#' @rdname expr_data
setMethod("expr_data", signature("CEMiTool"),
         function(cem_obj, filtered=T){
            if (filtered) {
                return(cem_obj@expression[cem_obj@selected_genes,])
            } else {
                return(cem_obj@expression)
            }
         })


#' @rdname expr_data
#' @export
setGeneric("expr_data<-", function(cem_obj, value) {
            standardGeneric("expr_data<-")
          })

#' @rdname expr_data
setReplaceMethod("expr_data", signature("CEMiTool"),
         function(cem_obj, value){
            cem_obj@expression <- value
            cem_obj@selected_genes <- rownames(cem_obj@expression)
            return(cem_obj)
         })



#' Retrieve and set mod_colors attribute
#'
#' @param cem_obj Object of class \code{CEMiTool}
#' @param value a character vector containing colors for each module
#'              the names should match with module names
#' 
#' @return 
#'
#' @rdname mod_colors
#' @export
setGeneric("mod_colors", function(cem_obj) {
            standardGeneric("mod_colors")
          })

#' @rdname mod_colors
setMethod("mod_colors", signature("CEMiTool"),
         function(cem_obj){
            mod_names <- unique(cem_obj@module$modules)
            nmod <- length(mod_names)
            cols <- cem_obj@mod_colors
            if(nmod != 0) { 
                if(length(cem_obj@mod_colors) == 0){
                    if(nmod <= 16) {
                        cols <- rainbow(16, s = 1, v = 0.7)[1:nmod]
                    } else {
                        cols <- rep(rainbow(16, s = 1, v = 0.7), ceiling(nmod/16))[1:nmod]
                    }
                    names(cols) <- mod_names
                } else {
                    if(!all(sort(names(cem_obj@mod_colors)) == sort(mod_names))){
                        warning("mod_colors does not match with modules!")
                    }
                }
            }
            return(cols)
         } )

#' @rdname mod_colors
#' @export
setGeneric("mod_colors<-", function(cem_obj, value) {
            standardGeneric("mod_colors<-")
          })

#' @rdname mod_colors
setReplaceMethod("mod_colors", signature("CEMiTool"),
         function(cem_obj, value){
            cem_obj@mod_colors <- value
            return(cem_obj)
         } )


#' Retrive or set the sample_annotation attribute
#'
#' @param cem_obj Object of class \code{CEMiTool}
#' @param value a data.frame containing the sample annotation, 
#'              should have at least two columns containing the Class 
#'              and the Sample Name that should match with samples in 
#'              expression
#'
#' @rdname sample_annotation
#' @export
setGeneric("sample_annotation", function(cem_obj) {
            standardGeneric("sample_annotation")
          })

#' @rdname sample_annotation
setMethod("sample_annotation", signature("CEMiTool"),
         function(cem_obj){
            return(cem_obj@sample_annotation)
         } )

#' @rdname sample_annotation
#' @export
setGeneric("sample_annotation<-", function(cem_obj, value) {
            standardGeneric("sample_annotation<-")
          })

#' @rdname sample_annotation
setReplaceMethod("sample_annotation", signature("CEMiTool"),
         function(cem_obj, value){
            if(!cem_obj@sample_name_column %in% colnames(value)){
                stop("Please supply a data.frame with a column named ", 
                     cem_obj@sample_name_column, 
                     " or change the slot sample_name_column.")
            }
            if(!cem_obj@class_column %in% colnames(value)){
                stop("Please supply a data.frame with a column named ", 
                     cem_obj@class_column, 
                     " or change the slot class_column.")
            }
            cem_obj@sample_annotation <- value
            return(cem_obj)
         } )


#' Full gene co-expression analysis
#'
#' Defines co-expression modules and functionally characterizes
#' each one of them.
#'
#' @param exprs Gene expression \code{data.frame}.
#' @param annot Sample annotation \code{data.frame}.
#' @param gmt a list from function prepare_gmt containing the gene sets.
#' @param interactions a data.frame containing two columns with gene names.
#' @param filter logical. If TRUE, will filter expression data.
#' @param filter_pval P-value threshold for filtering.
#' @param cor_method A character string indicating which correlation coefficient is
#'        to be computed. One of \code{"pearson"} or \code{"spearman"}.
#'        Default \code{"pearson"}.
#' @param sample_name_column A character string indicating the sample column 
#'        name of the annotation table.
#' @param class_column A character string indicating the class column name of the 
#'        annotation table.
#' @param merge_similar Logical. If \code{TRUE}, merge similar modules.
#' @param split_modules Logical. If \code{TRUE}, splits modules by correlation sign.
#' @param ora_pval P-value for overrepresentation analysis. Default \code{0.05}.
#' @param min_ngen Minimal number of genes per submodule. Default \code{30}.
#' @param diss_thresh Module merging correlation threshold for eigengene similarity.
#'        Default \code{0.8}.
#' @param plot Logical. If \code{TRUE}, plots all figures.
#' @param directed Logical. If \code{TRUE}, the igraph objects in interactions slot will be directed.
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
                     interactions,
                     filter=TRUE,
                     filter_pval=0.1,
                     cor_method=c('pearson', 'spearman'),
                     sample_name_column="SampleName",
                     class_column="Class",
                     merge_similar=TRUE,
                     split_modules=FALSE,
                     ora_pval=0.05,
                     min_ngen=30,
                     diss_thresh=0.8,
                     plot=FALSE,
                     directed=FALSE,
                     verbose=FALSE
                     )
{
    if (missing(exprs)) {
        stop('Must provide, at least, expression data')
    }
 
    # initialize CEMiTool object to hold data, analysis results and plots
    results <- new('CEMiTool', expression=exprs, 
                   sample_name_column=sample_name_column, 
                   class_column=class_column)
    
    if (filter) {
        results <- filter_expr(results, filter_pval)
        if (length(results@selected_genes) >= 0) {
            stop('Stopping analysis')
        }
    }

    results <- find_modules(results,
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

    if (!missing(interactions)){
        results <- include_interactions(results, int_df=interactions, 
                                        directed=directed)
    }

    # if user provides annot file
    if (!missing(annot)) {
        sample_annotation(results) <- annot
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

        if (!is.null(results@interactions)) {
            results <- plot_interactions(results)
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
#' @param object Object of class CEMiTool
#'
#' @export
setMethod('show', signature(object='CEMiTool'),
          function(object) {
              cat("CEMiTool Object\n")
              cat("- Number of modules:", nmodules(object), "\n")
              cat("- Modules: ")
              if(nrow(object@module) == 0){
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


