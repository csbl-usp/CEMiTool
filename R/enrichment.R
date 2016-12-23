#' Module Overrepresentation Analysis
#'
#' Perfoms overrepresentation analysis for each co-expression module found
#'
#' @param test
#'
#' @return None
#'
#' @examples
#' mod_ora(test)
#'
#' @export
mod_ora <- function(test){}


#' Module Gene Set Enrichment Analysis 
#'
#' Perfoms gene set enrichment analysis for each co-expression module found
#'
#' @param exprs 
#' @param gene_module
#'
#' @return 
#'
#' @examples
#' mod_gsea(test)
#'
#' @export
mod_gsea <- function(exprs, gene_module, annot, sample_col=1,
                     class_col=2)
{
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
        # samples of class == class_group
        class_samples <- annot[annot[, class_col]==class_group, sample_col]
        
        # genes ranked by mean
        genes_ranked <- apply(z_exprs[, class_samples], 1, mean)
        genes_ranked <- sort(genes_ranked, decreasing=TRUE)
        #register(SerialParam())
        gsea_results <- fgsea(pathways=gene_sets,
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



