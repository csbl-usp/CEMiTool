#' @import data.table
#' @import fgsea
NULL

#' Reads a GMT file
#'
#' @param fname GMT file name. 
#'
#' @return a list containing genes and description of each pathway
#'
#' @export
read_gmt <- function(fname){
    res <- list(genes=list(), desc=list())
    gmt <- file(fname)
    gmt_lines <- readLines(gmt)
    close(gmt)
    gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, split="\t")))
    gmt_names <- sapply(gmt_list, '[', 1)
    gmt_desc <- lapply(gmt_list, '[', 2)
    gmt_genes <- lapply(gmt_list, function(x){x[3:length(x)]})
    names(gmt_desc) <- names(gmt_genes) <- gmt_names
    res <- list()
    res[["term2gene"]] <- do.call(rbind, lapply(names(gmt_genes),
                                                function(n) cbind.data.frame(Term=n, Gene=gmt_genes[[n]], stringsAsFactors=FALSE)))
    res[["term2name"]] <- do.call(rbind, lapply(names(gmt_desc),
                                                function(n) cbind.data.frame(Term=n, Name=gmt_desc[[n]], stringsAsFactors=FALSE)))
    return(res)
}



# Performs Over Representation Analysis for a list of genes and a GMT
#
# @keywords internal
#
# @param topgenes a vector of genes
# @param gmt.list a gmt from prepare.gmt function
# @param allgenes a vector containing all genes to be considered as universe
#
# @return a data.frame containing the results
#
#
ora <- function(mod_name, gmt_list, allgenes, mods){
    if(missing(allgenes)) {
        message("Using all genes in GMT file as universe.")
        allgenes <- unique(gmt_list[["term2gene"]][, "Gene"])
    }
    topgenes <- mods[[mod_name]]
    enriched <- clusterProfiler::enricher(gene = topgenes,
                                          pvalueCutoff = 1,
                                          qvalueCutoff = 1,
                                          universe = allgenes,
                                          TERM2GENE = gmt_list[['term2gene']],
                                          TERM2NAME = gmt_list[['term2name']])
    if (!is.null(enriched) && !is.logical(enriched)) {
        result <- enriched@result
    } else {
        if(mod_name != "Not.Correlated"){
            warning(paste("Enrichment for module", mod_name, "is NULL"))    
        }
        result <- data.frame(Module=character(), ID=character(),
                             Description=character(),
                             GeneRation=numeric(), BgRatio=numeric(),
                             pvalue=numeric(), p.adjust=numeric(),
                             qvalue=numeric(), geneID=character(),
                             Count=numeric(), stringsAsFactors=FALSE)
    }
    return(result)
}


#' Module Overrepresentation Analysis
#'
#' Perfoms overrepresentation analysis for each co-expression module found.
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param gmt_in Output of CEMiTool::read_gmt function.
#' @param verbose logical. Report analysis steps.
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool}
#'
#' @examples
#' gmt <- read_gmt(system.file('extdata', 'pathways.gmt',
#'                    package='CEMiTool'))
#' mod_ora(cem, gmt)
#'
#' @rdname mod_ora
#' @export
setGeneric('mod_ora', function(cem, ...) {
    standardGeneric('mod_ora')
})

#' @rdname mod_ora
setMethod('mod_ora', signature(cem='CEMiTool'),
          function(cem, gmt_in, verbose=FALSE) {
              if (verbose) {
                  message('Running ORA')
              }
              message("Using all genes in GMT file as universe.")
              allgenes <- unique(gmt_in[["term2gene"]][, "Gene"])
              mods <- split(cem@module[, "genes"], cem@module[, "modules"])
              res_list <- lapply(names(mods), ora, gmt_in, allgenes, mods)
              if (all(lapply(res_list, nrow) == 0)){
                  message("Enrichment is NULL. Either your gmt file isn't good or your modules really aren't enriched for any of the pathways in the gmt file")
              }
              names(res_list) <- names(mods)
              for(n in names(res_list)) {
                  if (nrow(res_list[[n]]) > 0){
                      res_list[[n]] <- cbind.data.frame(Module=n, res_list[[n]], stringsAsFactors=FALSE)
                  }
              }
              res <- do.call(rbind, res_list)
              rownames(res) <- NULL
              cem@ora <- res
              return(cem)
          }
)

#' Module Gene Set Enrichment Analysis 
#'
#' Perfoms gene set enrichment analysis for each co-expression module found.
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param verbose logical. Report analysis steps.
#' @param ... Optional parameters.
#' 
#' @return GSEA results.
#'
#' @examples
#' mod_gsea(cem)
#'
#' @rdname mod_gsea
#' @export
setGeneric('mod_gsea', function(cem, ...) {
    standardGeneric('mod_gsea')
})

#' @rdname mod_gsea
setMethod('mod_gsea', signature(cem='CEMiTool'),
          function(cem, verbose=FALSE) {
              if (nrow(cem@sample_annotation)==0) {
                  warning('Looks like your sample_annotation slot is empty. Cannot proceed with gene set enrichment analysis.')
                  return(cem)
              }
              
              if (verbose) {
                  message('Running GSEA')
              }
              
              # creates gene sets from modules
              modules <- unique(cem@module[, 'modules'])
              gene_sets <- lapply(modules, function(mod){
                  return(cem@module[cem@module[, 'modules']==mod, 'genes'])
              })
              names(gene_sets) <- modules
              
              annot <- cem@sample_annotation
              class_col <- cem@class_column
              sample_col <- cem@sample_name_column
              classes <- unique(annot[, class_col])
              
              
              # Check if expression samples are all in sample_annotation
              expr_samples <- names(expr_data(cem))
              annot_samples <- sample_annotation(cem)[, sample_col]
              if(!all(expr_samples %in% annot_samples)){
                  stop("Sample annotation file does not contain all samples in expression file. Please input new sample annotation file using function sample_annotation()")
              }
              if(length(expr_samples) < length(annot_samples)){
                  warning("Expression file has less samples than annotation file. Cutting annotation file.")
                  annot <- annot[annot[,sample_col] %in% expr_samples,]
                  annot_samples <- annot[, sample_col]
              }
              
              # expression to z-score
              z_expr <- data.frame(t(scale(t(expr_data(cem, filtered=FALSE)), 
                                           center=TRUE, 
                                           scale=TRUE)),
                                   check.names = FALSE,
                                   stringsAsFactors=FALSE)
              
              # calculates enrichment for each module for each class in annot
              
              gsea_list <- lapply(classes, function(class_group){
                  
                  if (verbose) {
                      message(paste0('Calculating modules enrichment analysis for class ',
                                     class_group))
                  }
                  # samples of class == class_group
                  class_samples <- annot[annot[, class_col]==class_group, sample_col]
                  
                  # genes ranked by mean
                  genes_ranked <- apply(z_expr[, class_samples, drop=FALSE], 1, mean)
                  genes_ranked <- sort(genes_ranked, decreasing=TRUE)
                  
                  # BiocParallel setting up
                  BiocParallel::register(BiocParallel::SerialParam())
                  
                  gsea_results <- fgsea::fgsea(pathways=gene_sets,
                                               stats=genes_ranked,
                                               minSize=15,
                                               maxSize=500,
                                               nperm=10000,
                                               nproc=0)
                  setDF(gsea_results)
                  gsea_results[, 'leadingEdge'] <- unlist(lapply(gsea_results[, 'leadingEdge'],
                                                                 function(ledges){
                                                                     ledges <- paste(ledges, collapse=",")
                                                                 }))
                  columns <- colnames(gsea_results) 
                  colnames(gsea_results) <- c(columns[1], paste0(columns[-1], "_", class_group))
                  return(gsea_results)
              })
              # merging all classes gsea results into one data.frame
              all_classes_df <- Reduce(function(x,y) {
                  merge(x,y, all=TRUE, by='pathway')
              }, gsea_list) 
              
              # separating ES / NES / pval
              patterns <- list('es'='^ES_','nes'='^NES_', 'pval'='^pval_')
              out_gsea <- lapply(patterns, function(pattern) {
                  desired_stat <- all_classes_df[, c('pathway',
                                                     grep(pattern, colnames(all_classes_df),value=TRUE))]
                  colnames(desired_stat) <- gsub(pattern, '', colnames(desired_stat))
                  return(desired_stat)
              })
              
              names(out_gsea) <- names(patterns)
              cem@enrichment <- out_gsea 
              return(cem)
          })
