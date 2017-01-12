#' @import data.table

# Reads a GMT file
#
# @keywords internal
#
# @param fname GMT file name 
#
# @return a list containing genes and description of each pathway 
#
#
#
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
    return(list(genes=gmt_genes, desc=gmt_desc))
}


# Transforms a GMT list into inputs for enricher
#
# @keywords internal
#
# @param gmt a GMT list from read.gmt function
#
# @return a list containing term2gene and term2name
#
#
prepare_gmt <- function(gmt){
    res <- list()
    res[["term2gene"]] <- do.call(rbind, lapply(names(gmt[["genes"]]),
                            function(n) cbind(Term=n, Gene=gmt[["genes"]][[n]])))
    res[["term2name"]] <- do.call(rbind, lapply(names(gmt[["desc"]]),
                            function(n) cbind(Term=n, Name=gmt[["desc"]][[n]])))
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
ora <- function(topgenes, gmt_list, allgenes){
    if(missing(allgenes)) {
        message("Using all genes in GMT file as universe.")
        allgenes <- unique(gmt_list[["term2gene"]][, "Gene"])
    }
    enriched <- clusterProfiler::enricher(gene = topgenes,
                         pvalueCutoff = 1,
                         qvalueCutoff = 1,
                         universe = allgenes,
                         TERM2GENE = gmt_list[['term2gene']],
                         TERM2NAME = gmt_list[['term2name']])
    result <- enriched@result
    return(result)
}


#' Module Overrepresentation Analysis
#'
#' Perfoms overrepresentation analysis for each co-expression module found.
#'
#' @param gene_module Two column \code{data.frame}. First column with
#'        gene identifiers and second column with module information.
#' @param gmt GMT file name.
#' @param verbose logical. Report analysis steps.
#'
#' @return None
#'
#' @examples
#' mod_ora(gene_module, gmt)
#'
#' @export
mod_ora <- function(gene_module, gmt_in, verbose=FALSE)
{
    if (verbose) {
        message('Running ORA')
    }
    message("Using all genes in GMT file as universe.")
    allgenes <- unique(gmt_in[["term2gene"]][, "Gene"])
    mods <- split(gene_module[, "genes"], gene_module[, "modules"])
    res_list <- lapply(mods, ora, gmt_in, allgenes)
    names(res_list) <- names(mods)
    for(n in names(res_list)) {
        res_list[[n]] <- cbind.data.frame(Module=n, res_list[[n]])
    }
    res <- do.call(rbind, res_list)
    rownames(res) <- NULL
    return(res)
}


#' Module Gene Set Enrichment Analysis 
#'
#' Perfoms gene set enrichment analysis for each co-expression module found.
#'
#' @param exprs Gene expression \code{data.frame}.
#' @param gene_module Two column \code{data.frame}. First column with
#'        gene identifiers and second column with module information.
#' 
#' @return GSEA results.
#'
#' @examples
#' mod_gsea(test)
#'
#' @export
mod_gsea <- function(exprs, gene_module, annot, sample_col=1,
                     class_col=2, verbose=F)
{
    if (verbose) {
        message('Running GSEA')
    }
    
    # creates gene sets from modules
    modules <- unique(gene_module[, 'modules'])
    gene_sets <- lapply(modules, function(mod){
        return(gene_module[gene_module[, 'modules']==mod, 'genes'])
    })
    names(gene_sets) <- modules
    
    # expression to z-score
    z_exprs <- data.frame(t(scale(t(exprs), 
                                  center=TRUE, 
                                  scale=TRUE)))

    # calculates enrichment for each module for each class in annot
    classes <- unique(annot[, class_col])
    gsea_list <- lapply(classes, function(class_group){
        
        if (verbose) {
            message(paste0('Calculating modules enrichment analysis for class ',
                           class_group))
        }
        # samples of class == class_group
        class_samples <- annot[annot[, class_col]==class_group, sample_col]
        
        # genes ranked by mean
        genes_ranked <- apply(z_exprs[, class_samples], 1, mean)
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
                                 grep(pattern, colnames(all_classes_df),value=T))]
       colnames(desired_stat) <- gsub(pattern, '', colnames(desired_stat))
       return(desired_stat)
    })
    
    names(out_gsea) <- names(patterns)
    
    return(out_gsea)
}



