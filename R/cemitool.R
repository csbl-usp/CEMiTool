#' @importFrom grDevices rainbow
#' @importFrom utils write.table
#' @importFrom methods new 'slot<-' show

setOldClass('gg')
setOldClass('ggplot')
setOldClass('gtable')

#' An S4 class to represent the CEMiTool analysis.
#'
#' @slot expression Gene expression \code{data.frame}.
#' @slot sample_annotation Sample annotation \code{data.frame}.
#' @slot fit_indices \code{data.frame} containing scale-free model fit,
#'   soft-threshold and network parameters.
#' @slot selected_genes Character \code{vector} containing the names of genes
#'   selected for analysis
#' @slot module Genes in modules information \code{data.frame}.
#' @slot enrichment \code{list} with modules enrichment results for sample classes.
#' @slot ora Over-representation analysis results \code{data.frame}.
#' @slot interactions \code{list} containing gene interactions present in
#'   modules.
#' @slot interaction_plot list of ggplot graphs with module gene interactions.
#' @slot profile_plot list of ggplot graphs with gene expression profile per module.
#' @slot enrichment_plot ggplot graph for enrichment analysis results.
#' @slot beta_r2_plot ggplot graph with scale-free topology fit results for each
#'   soft-threshold.
#' @slot mean_k_plot ggplot graph with mean network connectivity.
#' @slot barplot_ora list of ggplot graphs with over-representation analysis results per module.
#' @slot sample_tree_plot gtable containing sample dendrogram with class labels and clinical data
#'      (if available in sample_annotation(cem)).
#' @slot mean_var_plot Mean x variance scatterplot.
#' @slot hist_plot Expression histogram.
#' @slot qq_plot Quantile-quantile plot.
#' @slot sample_name_column character string containing the name of the column with sample names
#'   in the annotation file.
#' @slot class_column character string containing the name of the column with class names
#'   in the annotation file.
#' @slot mod_colors character \code{vector} containing colors associated with each network module.
#' @slot parameters \code{list} containing analysis parameters.
#' @slot adjacency \code{matrix} containing gene adjacency values based on correlation
#'
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Initialize CEMiTool object with expression
#' cem <- new("CEMiTool", expression=expr0)
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
                                enrichment_plot='list',
                                beta_r2_plot='list',
                                mean_k_plot='list',
                                barplot_ora='list',
                                sample_tree_plot='gtable',
                                mean_var_plot='gg',
                                hist_plot='gg',
                                qq_plot='gg',
                                sample_name_column='vector',
                                class_column='vector',
                                mod_colors='character',
                                parameters='list',
                                input_params='list',
                                calls='list',
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

#' Create a CEMiTool object
#'
#' @param expr Object of class \code{data.frame} with gene
#'        expression data
#' @param sample_annot Object of \code{data.frame} containing the sample
#'        annotation. It should have at least two columns containing group Class
#'        and the Sample Name that should match with samples in
#'        expression file
#' @param sample_name_column A string specifying the column to be used as
#'        sample identification. Default: "SampleName".
#' @param class_column A string specifying the column to be used as
#'        a grouping factor for samples. Default: "Class"
#'
#' @return Object of class \code{CEMiTool}
#' @examples
#' # Create new CEMiTool object
#' cem <- new_cem()
#' # Create new CEMiTool object with expression and sample_annotation data
#' data(expr0)
#' data(sample_annot)
#' cem <- new_cem(expr0, sample_annot, "SampleName", "Class")
#' # Equivalent to a call to new()
#' cem2 <- new("CEMiTool", expression=expr0, sample_annotation=sample_annot)
#' identical(cem, cem2)
#' @export
new_cem <- function(expr=data.frame(), sample_annot=data.frame(),
        sample_name_column="SampleName", class_column="Class"){

    if(nrow(sample_annot) > 0){
        if(!sample_name_column %in% names(sample_annot)){
            stop("Please supply a data.frame with a column named ",
                sample_name_column,
                " or change the sample_name_column argument.")
            }
        if(!class_column %in% names(sample_annot)){
            stop("Please supply a data.frame with a column named ",
                class_column,
                " or change the class_column argument.")
        }
    }
    cem <- new("CEMiTool", expression=expr, sample_annotation=sample_annot,
        sample_name_column=sample_name_column, class_column=class_column)
    #cem <- get_args(cem, vars=mget(ls()))
    return(cem)
}

#' Retrieve and set expression attribute
#'
#' @param cem Object of class \code{CEMiTool}
#' @param value Object of class \code{data.frame} with gene
#'        expression data
#' @param filtered logical. If TRUE retrieves filtered expression data
#' @param ... Optional parameters.
#'
#' @return Object of class \code{data.frame} with gene expression data
#' @examples
#' # Initialize an empty CEMiTool object
#' cem <- new_cem()
#' # Get example expression data
#' data(expr0)
#' # Add expression file to CEMiTool object
#' expr_data(cem) <- expr0
#' # Check expression file
#' head(expr_data(cem))
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
        #cem <- get_args(cem, vars=mget(ls()))
        cem@expression <- value
        cem@selected_genes <- rownames(cem@expression)
        return(cem)
    })



#' Retrieve and set mod_colors attribute
#'
#' @param cem Object of class \code{CEMiTool}
#' @param value a character vector containing colors for each module.
#'Names should match with module names
#'
#' @return A vector with color names.
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # See module colors
#' mod_colors(cem)
#' @rdname mod_colors
#' @export
setGeneric("mod_colors", function(cem) {
    standardGeneric("mod_colors")
})

#' @rdname mod_colors
setMethod("mod_colors", signature("CEMiTool"),
    function(cem){
       mod_names <- mod_names(cem)
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
               if(is.null(names(cem@mod_colors))){
                   warning("mod_colors should be a character vector with names corresponding to the modules")
               } else if(!all(sort(names(cem@mod_colors)) == sort(mod_names))){
                   warning("mod_colors names do not match with modules!")
               }
           }
       }
       return(cols)
    })

#' @rdname mod_colors
#' @export
setGeneric("mod_colors<-", function(cem, value) {
    standardGeneric("mod_colors<-")
})

#' @rdname mod_colors
setReplaceMethod("mod_colors", signature("CEMiTool", "character"),
    function(cem, value){
        mod_names <- mod_names(cem)
        cem@mod_colors <- value
        if(is.null(names(cem@mod_colors))){
            stop("mod_colors should be a character vector with names corresponding to the modules")
        } else if(!all(sort(names(cem@mod_colors)) == sort(mod_names))){
            stop("mod_colors names do not match with modules!")
        }
        return(cem)
    })


#' Retrive or set the sample_annotation attribute
#'
#' @param cem Object of class \code{CEMiTool}
#' @param sample_name_column A string containing the name of a column which should be used
#' as a unique identifier for samples in the file. Only used when assigning a sample annotation
#' data.frame. Default: "SampleName".
#' @param class_column A string containing the name of a column which should be used
#' to identify different sample groups. Only used when assigning a sample annotation
#' data.frame. Default: "Class"
#' @param value A data.frame containing the sample annotation,
#' should have at least two columns containing the Class
#' and the Sample Name that should match with samples in
#' expression
#'
#' @return A data.frame containing characteristics of each sample.
#'
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Get example sample_annotation data
#' data(sample_annot)
#' # Initialize CEMiTool object with expression
#' cem <- new_cem(expr0)
#' # Add sample annotation file to CEMiTool object
#' sample_annotation(cem,
#'     sample_name_column="SampleName",
#'     class_column="Class") <- sample_annot
#' # Check annotation
#' head(sample_annotation(cem))
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
    })

#' @rdname sample_annotation
#' @export
setGeneric("sample_annotation<-", function(cem, sample_name_column="SampleName", class_column="Class", value) {
    standardGeneric("sample_annotation<-")
})

#' @rdname sample_annotation
setReplaceMethod("sample_annotation", signature("CEMiTool"),
    function(cem, sample_name_column="SampleName", class_column="Class", value){
       if(!sample_name_column %in% colnames(value)){
           stop("Please supply a data.frame with a column named ",
                sample_name_column,
                " or change the sample_name_column argument.")
       }
       if(!class_column %in% colnames(value)){
           stop("Please supply a data.frame with a column named ",
                class_column,
                " or change the class_column argument.")
       }
       if(min(table(value[, class_column])) == 1){
           warning("There is at least one class with only 1 sample in it. Results may be suboptimal.")
       }
       cem@sample_annotation <- value
       cem@sample_name_column <- sample_name_column
       cem@class_column <- class_column
       return(cem)
    })


#' Full gene co-expression analysis
#'
#' Defines co-expression modules and runs several different analyses.
#'
#' @param expr Gene expression \code{data.frame}.
#' @param annot Sample annotation \code{data.frame}.
#' @param gmt A data.frame containing two columns, one with pathways and one with genes
#' @param interactions A data.frame containing two columns with gene names.
#' @param filter Logical. If TRUE, will filter expression data.
#' @param filter_pval P-value threshold for filtering.Default \code{0.1}.
#' @param apply_vst Logical. If TRUE, will apply Variance Stabilizing Transform before filtering genes.
#'           Currently ignored if parameter \code{filter} is FALSE.
#' @param n_genes Number of genes left after filtering.
#' @param eps A value for accepted R-squared interval between subsequent beta values. Default is 0.1. 
#' @param cor_method A character string indicating which correlation coefficient is
#'        to be computed. One of \code{"pearson"} or \code{"spearman"}.
#'        Default is \code{"pearson"}.
#' @param cor_function A character string indicating the correlation function to be used. Supported functions are
#'           currently 'cor' and 'bicor'. Default is \code{"cor"}
#' @param network_type A character string indicating if network type should be computed
#'        as \code{"signed"} or \code{"unsigned"}. Default is \code{"unsigned"}
#' @param tom_type  A character string indicating if the TOM type should be computed
#'        as \code{"signed"} or \code{"unsigned"}. Default is \code{"signed"}
#' @param set_beta A value to override the automatically selected beta value. Default is NULL.
#' @param force_beta Whether or not to automatically force a beta value based on number of samples. Default is FALSE.
#' @param sample_name_column A character string indicating the sample column
#'        name of the annotation table.
#' @param class_column A character string indicating the class column name of the
#'        annotation table.
#' @param merge_similar Logical. If \code{TRUE}, merge similar modules.
#' @param rank_method Character string indicating how to rank genes. Either "mean" 
#'        (the default) or "median".
#' @param ora_pval P-value for overrepresentation analysis. Default \code{0.05}.
#' @param gsea_scale If TRUE, apply z-score transformation for GSEA analysis. Default is \code{TRUE}
#' @param gsea_min_size Minimum size of gene sets for GSEA analysis. Default is \code{15}
#' @param gsea_max_size Maximum size of gene sets for GSEA analysis. Default is \code{500}
#' @param min_ngen Minimal number of genes per submodule. Default \code{30}.
#' @param diss_thresh Module merging correlation threshold for eigengene similarity.
#'        Default \code{0.8}.
#' @param plot Logical. If \code{TRUE}, plots all figures.
#' @param order_by_class Logical. If \code{TRUE}, samples in profile plot are ordered by the groups  
#'           defined by the class_column slot in the sample annotation file. Ignored if there is no 
#'           sample_annotation file. Default \code{TRUE}.
#' @param center_func Character string indicating the centrality measure to show in
#'        the plot. Either 'mean' (the default) or 'median'.
#' @param directed Logical. If \code{TRUE}, the igraph objects in interactions slot will be directed.
#' @param verbose Logical. If \code{TRUE}, reports analysis steps.
#'
#' @return Object of class \code{CEMiTool}
#'
#' @examples
#' \dontrun{
#' # Get example expression data
#' data(expr0)
#' # Run CEMiTool analyses
#' cem <- cemitool(expr=expr0)
#' # Run CEMiTool applying Variance Stabilizing Transformation to data
#' cem <- cemitool(expr=expr0, apply_vst=TRUE)
#' # Run CEMiTool with additional processing messages
#' cem <- cemitool(expr=expr0, verbose=TRUE)
#'
#' # Run full CEMiTool analysis
#' ## Get example sample annotation data
#' data(sample_annot)
#' ## Read example pathways file
#' gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
#' gmt_in <- read_gmt(gmt_fname)
#' ## Get example interactions file
#' int_df <- read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool"))
#' ## Run CEMiTool
#' cem <- cemitool(expr=expr0, annot=sample_annot, gmt=gmt_in,
#'     interactions=int_df, verbose=TRUE, plot=TRUE)
#'
#' # Create report as html file
#' generate_report(cem, directory = "./Report", output_format="html_document")
#'
#' # Write analysis results into files
#' write_files(cem, directory="./Tables", force=TRUE)
#' 
#' # Save all plots
#' save_plots(cem, "all", directory="./Plots")
#' }
#' 
#' 
#' @export
cemitool <- function(expr,
                     annot,
                     gmt,
                     interactions,
                     filter=TRUE,
                     filter_pval=0.1,
                     apply_vst=FALSE,
                     n_genes,
                     eps=0.1,
                     cor_method=c('pearson', 'spearman'),
                     cor_function='cor',
                     network_type='unsigned',
                     tom_type='signed',
                     set_beta=NULL,
                     force_beta=FALSE, 
                     sample_name_column="SampleName",
                     class_column="Class",
                     merge_similar=TRUE,
                     rank_method="mean",
                     ora_pval=0.05,
                     gsea_scale=TRUE,
                     gsea_min_size=15,
                     gsea_max_size=500,
                     min_ngen=30,
                     diss_thresh=0.8,
                     plot=TRUE,
                     order_by_class=TRUE,
                     center_func="mean",
                     directed=FALSE,
                     verbose=FALSE
                     ){

    if (missing(expr)) {
        stop('Please provide expression data!')
    }

    # initialize CEMiTool object to hold data, analysis results and plots
    results <- new('CEMiTool', expression=expr,
                   sample_name_column=sample_name_column,
                   class_column=class_column)
    
    # keep input parameters
    #results <- get_args(cem=results, vars=mget(ls()))

    if (filter) {
        if(!missing(n_genes)){
            results <- filter_expr(results, n_genes=n_genes, apply_vst=apply_vst)
        } else {
            results <- filter_expr(results, filter_pval=filter_pval, apply_vst=apply_vst)
        }
        if (length(results@selected_genes) <= 0) {
            stop('Stopping analysis, no gene left for analysis. Maybe try to change the filter parameters.')
        }
    }

    # if user provides annot file
    if (!missing(annot)) {
        if(verbose){
            message("Including sample annotation ...")
        }
        sample_annotation(results, sample_name_column=sample_name_column, class_column=class_column) <- annot
    }

    if(plot){
        if(verbose){
            message("Plotting diagnostic plots ...")
            message("...Plotting mean and variance scatterplot ...")
            message("...Plotting expression histogram ...")
            message("...Plotting qq plot ...")
            message("...Plotting sample tree ...")
        }
        results <- plot_mean_var(results)
        results <- plot_hist(results)
        results <- plot_qq(results)
        results <- plot_sample_tree(results)
    }

    if(verbose){
        message("Finding modules ...")
    }
    results <- find_modules(results,
                            cor_method=match.arg(cor_method),
                            cor_function=cor_function,
                            eps=0.1,
                            min_ngen=min_ngen,
                            merge_similar=merge_similar,
                            diss_thresh=diss_thresh,
                            network_type=network_type,
                            tom_type=tom_type,
                            set_beta=set_beta,
                            force_beta=force_beta,
                            verbose=verbose)
    if(verbose){
        message("Plotting beta x R squared curve ...")
        message("Plotting mean connectivity curve ...")
    }

    results <- plot_beta_r2(results)
    results <- plot_mean_k(results)

    if(is.na(results@parameters$beta)){
        message("Unable to find parameter beta. Please check diagnostic plots with function diagnostic_report().")
        return(results)
    }

    if (!missing(interactions)){
        if(verbose){
            message("Including interactions ...")
        }
        interactions_data(results) <- interactions
    }


    if (!missing(annot)){
        if(verbose){
            message("Running Gene Set Enrichment Analysis ...")
        }
        #run mod_gsea
        results <- mod_gsea(results, gsea_scale=gsea_scale, rank_method=rank_method, 
                            gsea_min_size=gsea_min_size, gsea_max_size=gsea_max_size, verbose=verbose)
    }

    # if user provides .gmt file
    if (!missing(gmt)) {
        if(verbose){
            message("Running over representation analysis ...")
        }
        #run mod_ora
        results <- mod_ora(results, gmt=gmt, verbose=verbose)
    }

    # plots all desired charts
    if (plot) {
        if(verbose){
            message("Generating profile plots ...")
        }
        results <- plot_profile(results, order_by_class=order_by_class, center_func=center_func)

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

            results <- plot_ora(results, pv_cut=ora_pval)
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


    results@calls <- results@calls["cemitool", drop=FALSE]
    results@input_params <- results@input_params["cemitool", drop=FALSE]

    #results@input_params <- as.list(results@input_params$cemitool)
    #results@calls <- as.list(results@calls$cemitool)
    
    #cem@input_params <- list()
    #cem@calls <- list()
    #results <- get_args(cem=results, vars=mget(ls()))

    return(results)
}


#' Get the number of modules in a CEMiTool object
#'
#' @param cem Object of class \code{CEMiTool}
#'
#' @return number of modules
#'
#' @rdname nmodules
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Get the number of modules
#' nmodules(cem)
#'
#' @export
setGeneric('nmodules', function(cem) {
    standardGeneric('nmodules')
})

#' @rdname nmodules
setMethod('nmodules', signature(cem='CEMiTool'),
    function(cem) {
        n <- 0
        if(nrow(cem@module) > 0){
            n <- length(unique(cem@module$modules))
        } else {
            warning("Run cemitool function to get modules!")
        }
        return(n)
    })

#' Get module names in a CEMiTool object
#'
#' @param cem Object of class \code{CEMiTool}
#' @param include_NC Logical. Whether or not to include "Not.Correlated"
#' module. Defaults to \code{TRUE}.
#'
#' @return Module names
#'
#' @rdname mod_names
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Get module names
#' mod_names(cem)
#'
#' @export
setGeneric('mod_names', function(cem, include_NC=TRUE) {
    standardGeneric('mod_names')
})

#' @rdname mod_names
setMethod('mod_names', signature(cem='CEMiTool'),
          function(cem, include_NC=TRUE) {
              mods <- NULL
              if(nrow(cem@module) > 0){
                  mods <- names(sort(table(cem@module$modules), decreasing=TRUE))
                  if(!include_NC && ("Not.Correlated" %in% mods)){
                      mods <- mods[mods != "Not.Correlated"]
                  }
              } else {
                  warning("No modules in this CEMiTool object.")
              }
              return(mods)
          }
)


#' Get the module genes in a CEMiTool object
#'
#' @param cem Object of class \code{CEMiTool}
#' @param module A character string with the name of the module of which
#' genes are to be returned. Defaults to \code{NULL}, which returns the full
#' list of genes and modules.
#'
#' @return Object of class \code{data.frame} containing genes and their
#' respective module
#'
#' @rdname module_genes
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Get the module genes
#' module_genes(cem)
#' # Get genes for module M1
#' module_genes(cem, module="M1")
#' @export
setGeneric('module_genes', function(cem, module=NULL) {
    standardGeneric('module_genes')
})

#' @rdname module_genes
setMethod('module_genes', signature(cem='CEMiTool'),
    function(cem, module=NULL){
        #mod_names <- unique(cem@module[, "modules"])
        res <- NULL
        if(nrow(cem@module) > 0){
            res <- cem@module
        }else{
            message("No modules in this CEMiTool object.")
            return(res)
        }
        mod_names <- unique(cem@module[, "modules"])
        if(!is.null(module)){
            if(module %in% mod_names){
                res <- res[res$modules==module,]
            }else{
                stop("Undefined module!")
            }
        }
        return(res)
    }
)

#' Print a cemitool object
#'
#' @param object Object of class CEMiTool
#'
#' @return A CEMiTool object.
#'
#' @export
setMethod('show', signature(object='CEMiTool'),
    function(object) {
        cat("CEMiTool Object\n")
        cat("- Number of modules:", suppressWarnings(nmodules(object)), "\n")
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

#' Transform module genes list to a gmt file
#'
#' @keywords internal
#' 
#' @param cem 
#'
#' @return A .gmt file containing module genes in each row
#'
module_to_gmt <- function(cem, directory="./Tables"){
    if(nrow(cem@module) == 0){
        stop("No modules in CEMiTool object! Did you run find_modules()?")
    }else{
        gene_modules <- cem@module
        n_genes <- as.numeric(table(gene_modules[, "modules"]))
        n_genes <- n_genes[1:(length(n_genes)-1)]
        module_names <- as.character(unique(gene_modules[, "modules"]))
        module_names <- module_names[-which(module_names =="Not.Correlated")]
        module_names <- module_names[order(nchar(module_names), module_names)]
        
        gmt_df  <- as.data.frame(matrix("", ncol = max(n_genes), nrow = length(n_genes)), stringsAsFactors = FALSE)
        
        rownames(gmt_df) <- module_names
        
        for (i in 1:length(module_names)){
            mod <- module_names[i]
            selected <- gene_modules[gene_modules$modules == mod, "genes"]
            gmt_df[mod, 1:(length(selected))] <- selected
        }
        
        gmt_df <- as.data.frame(cbind(module_names, gmt_df))
        write.table(gmt_df, file.path(directory, "modules_genes.gmt"), sep="\t", col.names = FALSE)
    }
}


#' Save the CEMiTool object in files
#'
#' @param cem Object of class \code{CEMiTool}
#' @param directory a directory
#' @param force if the directory exists the execution will not stop
#' @param ... Optional parameters
#' @return A directory containing CEMiTool results in files.
#' @examples
#' \dontrun{
#' # Get example CEMiTool object
#' data(cem)
#' # Save CEMiTool results in files
#' write_files(cem, directory=".", force=TRUE)
#' }
#' 
#' @rdname write_files
#' @export
setGeneric('write_files', function(cem, ...) {
    standardGeneric('write_files')
})

#' @rdname write_files
setMethod('write_files', signature(cem='CEMiTool'),
    function(cem, directory="./Tables", force=FALSE) {
        if(dir.exists(directory)){
            if(!force){
                stop("Stopping analysis: ", directory, " already exists! Use force=TRUE to overwrite.")
            }
        } else {
            dir.create(directory, recursive=TRUE)
        }
        if(nrow(cem@module) > 0){
            write.table(cem@module, file.path(directory, "module.tsv"), sep="\t", row.names=FALSE)

            mean_summary <- mod_summary(cem, "mean")
            write.table(mean_summary, file.path(directory, "summary_mean.tsv"), sep="\t", row.names=FALSE)

            median_summary <- mod_summary(cem, "median")
            write.table(median_summary, file.path(directory, "summary_median.tsv"), sep="\t", row.names=FALSE)

            eg_summary <- mod_summary(cem, "eigengene")        
            write.table(eg_summary, file.path(directory, "summary_eigengene.tsv"), sep="\t", row.names=FALSE)
            
            module_to_gmt(cem, directory=directory)
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
            mod_ints <- lapply(names(cem@interactions), function(x){
                mod_int <- igraph::get.edgelist(cem@interactions[[x]])
                if(nrow(mod_int) > 0 ){
                    cbind(x, mod_int)
                }
            })
            int_df <- do.call("rbind", mod_ints)
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
