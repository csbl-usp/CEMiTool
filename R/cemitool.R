#' @importFrom grDevices rainbow
#' @importFrom utils head write.table
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
                                fit_indices='data.frame',
                                selected_genes='vector',
                                module='data.frame',
                                enrichment='list', # gene set enrichment analysis
                                ora='data.frame',
                                interactions='list',
                                interaction_plot='list',
                                profile_plot='list',
                                enrichment_plot='gg',
                                beta_r2_plot='gg',
                                mean_k_plot='gg',
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
#' @param cem Object of class \code{CEMiTool}
#' @param value Object of class \code{data.frame} with gene
#'        expression data
#' @param filtered logical. If TRUE retrieves filtered expression data
#' @param ... Optional parameters.
#'
#' @return Object of class \code{data.frame} with gene expression data
#'
#' @rdname expr_data
#' @export
setGeneric("expr_data", function(cem, ...) {
            standardGeneric("expr_data")
          })

#' @rdname expr_data
setMethod("expr_data", signature("CEMiTool"),
         function(cem, filtered=TRUE){
            if (filtered) {
                return(cem@expression[cem@selected_genes,])
            } else {
                return(cem@expression)
            }
         })


#' @rdname expr_data
#' @export
setGeneric("expr_data<-", function(cem, value) {
            standardGeneric("expr_data<-")
          })

#' @rdname expr_data
setReplaceMethod("expr_data", signature("CEMiTool"),
         function(cem, value){
            cem@expression <- value
            cem@selected_genes <- rownames(cem@expression)
            return(cem)
         })



#' Retrieve and set mod_colors attribute
#'
#' @param cem Object of class \code{CEMiTool}
#' @param value a character vector containing colors for each module
#'              the names should match with module names
#' 
#' @rdname mod_colors
#' @export
setGeneric("mod_colors", function(cem) {
            standardGeneric("mod_colors")
          })

#' @rdname mod_colors
setMethod("mod_colors", signature("CEMiTool"),
         function(cem){
            mod_names <- unique(cem@module$modules)
            nmod <- length(mod_names)
            cols <- cem@mod_colors
            if(nmod != 0) {
                if(length(cem@mod_colors) == 0){
                    if(nmod <= 16) {
                        cols <- rainbow(16, s = 1, v = 0.7)[1:nmod]
                    } else {
                        cols <- rep(rainbow(16, s = 1, v = 0.7), ceiling(nmod/16))[1:nmod]
                    }
                    names(cols) <- mod_names
                } else {
                    if(!all(sort(names(cem@mod_colors)) == sort(mod_names))){
                        warning("mod_colors does not match with modules!")
                    }
                }
            }
            return(cols)
         } )

#' @rdname mod_colors
#' @export
setGeneric("mod_colors<-", function(cem, value) {
            standardGeneric("mod_colors<-")
          })

#' @rdname mod_colors
setReplaceMethod("mod_colors", signature("CEMiTool"),
         function(cem, value){
            cem@mod_colors <- value
            return(cem)
         } )


#' Retrive or set the sample_annotation attribute
#'
#' @param cem Object of class \code{CEMiTool}
#' @param value a data.frame containing the sample annotation,
#'              should have at least two columns containing the Class
#'              and the Sample Name that should match with samples in
#'              expression
#'
#' @rdname sample_annotation
#' @export
setGeneric("sample_annotation", function(cem) {
            standardGeneric("sample_annotation")
          })

#' @rdname sample_annotation
setMethod("sample_annotation", signature("CEMiTool"),
         function(cem){
            return(cem@sample_annotation)
         } )

#' @rdname sample_annotation
#' @export
setGeneric("sample_annotation<-", function(cem, value) {
            standardGeneric("sample_annotation<-")
          })

#' @rdname sample_annotation
setReplaceMethod("sample_annotation", signature("CEMiTool"),
         function(cem, value){
            if(!cem@sample_name_column %in% colnames(value)){
                stop("Please supply a data.frame with a column named ",
                     cem@sample_name_column,
                     " or change the slot sample_name_column.")
            }
            if(!cem@class_column %in% colnames(value)){
                stop("Please supply a data.frame with a column named ",
                     cem@class_column,
                     " or change the slot class_column.")
            }
            if(min(table(value[, cem@class_column])) == 1){
                warning("There is at least one class with only 1 sample in it. Results may be suboptimal.")
            }
            cem@sample_annotation <- value
            return(cem)
         } )


#' Full gene co-expression analysis
#'
#' Defines co-expression modules and functionally characterizes
#' each one of them.
#'
#' @param expr Gene expression \code{data.frame}.
#' @param annot Sample annotation \code{data.frame}.
#' @param gmt A list from function read_gmt containing the gene sets.
#' @param interactions A data.frame containing two columns with gene names.
#' @param filter Logical. If TRUE, will filter expression data.
#' @param filter_pval P-value threshold for filtering.Default \code{0.1}.
#' @param n_genes Number of genes left after filtering.
#' @param cor_method A character string indicating which correlation coefficient is
#'        to be computed. One of \code{"pearson"} or \code{"spearman"}.
#'        Default is \code{"pearson"}.
#' @param network_type A character string indicating if network type should be computed 
#'        as \code{"signed"} or \code{"unsigned"}. Default is \code{"unsigned"}
#' @param tom_type  A character string indicating if the TOM type should be computed 
#'        as \code{"signed"} or \code{"unsigned"}. Default is \code{"signed"}
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
#' @return Object of class \code{CEMiTool}
#'
#' @examples
#' cemitool(expr=expr)
#'
#' @export
cemitool <- function(expr,
                     annot,
                     gmt,
                     interactions,
                     filter=TRUE,
                     filter_pval=0.1,
                     n_genes,
                     cor_method=c('pearson', 'spearman'),
                     network_type='unsigned',
                     tom_type='signed',
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
    if (missing(expr)) {
        stop('Must provide, at least, expression data')
    }

    # initialize CEMiTool object to hold data, analysis results and plots
    results <- new('CEMiTool', expression=expr,
                   sample_name_column=sample_name_column,
                   class_column=class_column)

    if (filter) {
        if(!missing(n_genes)){
            results <- filter_expr(results, n=n_genes)
        } else {
            results <- filter_expr(results, filter_pval)
        }
        if (length(results@selected_genes) <= 0) {
            stop('Stopping analysis, no gene left for analysis, try to change the filter parameters.')
        }
    }

    if(verbose){
        message("Finding modules ...")
    }
    results <- find_modules(results,
                            cor_method=match.arg(cor_method),
                            min_ngen=min_ngen,
                            merge_similar=merge_similar,
                            diss_thresh=diss_thresh,
                            network_type=network_type,
                            tom_type=tom_type,
                            verbose=verbose)

    # if user wants splitted modules
    if (split_modules) {
        if(verbose){
            message("Spliting modules ...")
        }
        results <- split_modules(results, min_ngen=min_ngen,
                                 verbose=verbose)
    }

    if (!missing(interactions)){
        if(verbose){
            message("Including interactions ...")
        }
        results <- include_interactions(results, int_df=interactions,
                                        directed=directed)
    }

    # if user provides annot file
    if (!missing(annot)) {
        if(verbose){
            message("Including sample annotation ...")
        }
        sample_annotation(results) <- annot
        if(verbose){
            message("Running Gene Set Enrichment Analysis ...")
        }
        #run mod_gsea
        results <- mod_gsea(results, verbose=verbose)
    }

    # if user provides .gmt file
    if (!missing(gmt)) {
        if(verbose){
            message("Running over representation analysis ...")
        }
        #run mod_ora
        results <- mod_ora(results, gmt_in=gmt, verbose=verbose)
    }

    # plots all desired charts
    if (plot) {
        if(verbose){
            message("Generating profile plots ...")
        }
        results <- plot_profile(results)

        if (length(results@enrichment) > 0) {
            if(verbose){
                message("Plotting Enrichment Scores ...")
            }
            results <- plot_gsea(results)
        }

        if (nrow(results@ora) > 0) {
            if(verbose){
                message("Plotting over representation analysis results ...")
            }

            results <- plot_ora(results)
        }

        if (length(results@interactions) > 0) {
            if(verbose){
                message("Plotting interaction network ...")
            }
            results <- plot_interactions(results)
        }
        
        if(verbose){
            message("Plotting beta x R squared curve ...")
            message("Plotting mean connectivity curve ...")
        }
      
        results <- plot_beta_r2(results)
        results <- plot_mean_k(results)
    }
    return(results)
}


#' Get the number of modules on a cemitool object
#'
#' @param cem Object of class \code{CEMiTool}
#'
#' @return number of modules
#'
#' @rdname nmodules
#' @export
setGeneric('nmodules', function(cem) {
    standardGeneric('nmodules')
})

#' @rdname nmodules
setMethod('nmodules', signature(cem='CEMiTool'),
          function(cem) {
              n <- 0
              if(!is.null(cem@module)){
                  n <- length(unique(cem@module$modules))
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
              cat("- Expression file: ")
              if(nrow(object@expression) == 0){
                  cat("null\n")
              } else {
                  cat("data.frame with", nrow(object@expression), "genes and", ncol(object@expression), "samples\n")
              }
              if(is.character(object@selected_genes)){
                  cat("- Selected data:", length(object@selected_genes), "genes selected\n")
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
              cat("- ORA barplot: ")
              if(length(object@barplot_ora)==0) {
                  cat("null\n")
              } else {
                  cat("ok\n")
              }
              cat("- Beta x R2 plot: ")
              if(!is.ggplot(object@beta_r2_plot)) {
                  cat("null\n")
              } else {
                  cat("ok\n")
              }
              cat("- Mean connectivity plot: ")
              if(!is.ggplot(object@mean_k_plot)) {
                  cat("null\n")
              } else {
                  cat("ok\n")
              }
          }
)




#' Save the CEMiTool object in files
#'
#' @param cem Object of class \code{CEMiTool}
#' @param directory a directory
#' @param force if the directory exists the execution will not stop
#' @param ... Optional parameters
#'
#' @examples 
#' write_files(cem, directory=".", force=TRUE)
#'
#' @rdname save_files
#' @export
setGeneric('write_files', function(cem, ...) {
    standardGeneric('write_files')
})

#' @rdname save_files
setMethod('write_files', signature(cem='CEMiTool'),
          function(cem, directory, force=FALSE) {
              if(dir.exists(directory)){
                  if(!force){
                      stop("Stopping analysis: ", directory, " already exists!")
                  } else {
                      warning(directory, " already exists!")
                  }
              } else {
                  dir.create(directory)
              }
              if(nrow(cem@module) > 0){
                  write.table(cem@module, file.path(directory, "module.tsv"), sep="\t", row.names=FALSE)
              }
              if(length(cem@selected_genes) > 0){
                  writeLines(cem@selected_genes, file.path(directory, "selected_genes.txt"))
              }
              if(length(cem@enrichment) > 0){
                  write.table(cem@enrichment$nes, file.path(directory, "enrichment_nes.tsv"), sep="\t", row.names=FALSE)
                  write.table(cem@enrichment$es, file.path(directory, "enrichment_es.tsv"), sep="\t", row.names=FALSE)
                  write.table(cem@enrichment$pval, file.path(directory, "enrichment_pval.tsv"), sep="\t", row.names=FALSE)
              }
              if(nrow(cem@ora) > 0){
                  write.table(cem@ora, file.path(directory, "ora.tsv"), sep="\t", row.names=FALSE)
              }
              if(length(cem@interactions) > 0){
                  int_df <- data.frame(Module=character(),
                                       Gene1=character(),
                                       Gene2=character())
                  for(n in names(cem@interactions)){
                      mod_int <- igraph::get.edgelist(cem@interactions[[n]])
                      if(nrow(mod_int) > 0 ){
                          int_df <- rbind.data.frame(int_df,
                                                 data.frame(Module=n,
                                                            igraph::get.edgelist(cem@interactions[[n]])
                                                            ))
                      }
                  }
                  colnames(int_df) <- c("Module", "Gene1", "Gene2")
                  write.table(int_df, file.path(directory, "interactions.tsv"), sep="\t", row.names=FALSE)
              }
              if(length(cem@parameters) > 0){
                  params <- cem@parameters
                  param_df <- data.frame(Parameter=names(params), Value=as.character(params))
                  write.table(param_df, file.path(directory, "parameters.tsv"), sep="\t", row.names=FALSE)
              }
          }
)
