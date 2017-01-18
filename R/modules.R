#' @importFrom grDevices pdf
#' @importFrom grDevices colorRampPalette
#' @import stats
#' @import data.table

.datatable.aware=TRUE
#' Co-expression modules definition
#'
#' Defines co-expression modules
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param cor_method A character string indicating which correlation coefficient
#'                   is to be computed.
#' @param min_ngen Minimal number of genes per submodule. Default \code{30}.
#' @param merge_similar Logical. If \code{TRUE}, merge similar modules.
#' @param diss_thresh Module merging correlation threshold for eigengene similarity.
#'        Default \code{0.8}.
#' @param verbose Logical. If \code{TRUE}, reports analysis steps.
#'
#' @return Object of class \code{CEMiTool} 
#'
#' @examples
#' data(exprs)
#' cem_obj <- new('CEMiTool', expression=exprs)
#' cem_obj <- find_modules(cem_obj)
#'
#' @rdname find_modules
#' @export
setGeneric('find_modules', function(cem_obj, ...) {
    standardGeneric('find_modules')
})

#' @rdname find_modules
#' @export
setMethod('find_modules', signature('CEMiTool'), 
          function(cem_obj, cor_method=c('pearson', 'spearman'),
                             min_ngen=30, merge_similar=T,
                             diss_thresh=0.8, verbose=F)
{

    if (is.null(exprs)){
        stop('Must provide expression data')
    }

    exprs_t <- t(cem_obj@expression)
    names(exprs_t) <- rownames(cem_obj@expression)
    rownames(exprs_t) <- colnames(cem_obj@expression)


    if (verbose) {
        message('Selecting Beta')
        verbosity <- 10
    } else {
        verbosity <- 0
    }

    # Define a range of soft-thresholding candidates
    powers <- c(c(1:10), seq(12, 20, 2))
    
    ## Automatic selection of soft-thresholding power beta ##
    beta <- WGCNA::pickSoftThreshold(exprs_t, powerVector=powers,
                              networkType='signed', moreNetworkConcepts=TRUE,
                              corOptions=list(use='p',
                                              method=match.arg(cor_method)),
                              verbose=verbosity)

    fit <- -sign(beta$fitIndices[, 3])*beta$fitIndices[, 2]
    k <- beta$fitIndices[, 5]
    At  <- powers[length(powers)] - powers[1]
    A   <- 0.0
    for(cont in 2:length(fit)) {
        A <- A + 0.5*(fit[cont] +
                      fit[cont-1])*(powers[cont] - powers[cont-1])
    }

    # Area under the curve/threshold
    phi <- A/At
    eps <- 0.1
    st <- c(NA,NA)

    # Selecting smallest beta value in Cauchy sequence range (CEMiTool beta)
    for(count in (1:(length(fit)-2))) {
        if(fit[count] >= 0.8) {
            d <- c(abs(fit[count] - fit[count+1]),
                   abs(fit[count] - fit[count+2]),
                   abs(fit[count+1] - fit[count+2]))
            if(max(d) < eps) {
                j <- which.max(k[count:count+2]) + count - 1
                st <- c(fit[j], powers[j])
                break
            }
        }
    }

    # Get network connectivities
    our_k <- NA
    if(!is.na(st[2])) {
        if(st[2] <= 10) {
            our_k <- beta$fitIndices[(st[2]) , 5]
        } else {
            line <- (st[2] + 10)/2
            our_k <- beta$fitIndices[line, 5]
        }
    }
    our_beta <- as.integer(st[2])
    our_r2 <- st[1]

    # Calculating adjacency matrix
    our_adj <- WGCNA::adjacency(exprs_t, power=our_beta, type='signed')

    # Calculating Topological Overlap Matrix
    our_tom <- WGCNA::TOMsimilarity(our_adj)

    # Determining TOM based distance measure
    our_diss <- 1 - our_tom

    # Clustering
    our_tree <- hclust(as.dist(our_diss), method = 'average')

    # Cutting tree to determine modules
    our_mods <- dynamicTreeCut::cutreeDynamic(dendro = our_tree, distM = our_diss,
                              deepSplit = 2,
                              pamRespectsDendro = FALSE,
                              minClusterSize = min_ngen)
    
    our_colors <- paste0('M',our_mods)
    
    #our_table <- table(our_mods)

    #our_colors <- WGCNA::labels2colors(our_mods)

    n_mods <- length(unique(our_colors))
    
    out <- data.frame(genes=rownames(exprs), modules=our_colors)

    # checks number of modules found
    if (n_mods <= 1) {
        stop('Could not specify the parameter Beta. No modules found.')
    }

    # if merge_similar=TRUE, merges similar modules
    if (merge_similar) {
        if (verbose) {
            message('Merging modules based on eigengene similarity')
        }
        # Calculates eigengenes
        me_list <- WGCNA::moduleEigengenes(exprs_t, colors=our_colors)
        me_eigen <- me_list$eigengenes

        # Calculates dissimilarity of module eigengenes
        me_diss <- 1 - cor(me_eigen)

        # Clustering module eigengenes
        me_tree <- hclust(as.dist(me_diss), method='average')

        # Setting cut height
        me_diss_thresh <- 1 - diss_thresh 

        # Merging modules
        merged_mods <-  WGCNA::mergeCloseModules(exprs_t, our_colors,
                                          cutHeight=diss_thresh)
        
        # The merged modules colors
        merged_mods <- factor(merged_mods$colors)

        levels(merged_mods) <- paste0('M', seq(1, length(unique(merged_mods)))) 

        out[, 'modules'] <- as.character(merged_mods)
    }
 
    cem_obj@module <- out
    return(cem_obj)
})


#' Co-expression module refinement 
#'
#' Refines modules by splitting them into submodules based on correlation sign.
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param min_ngen Minimal number of genes per submodule. Default \code{30}.
#' @param verbose Logical. If \code{TRUE}, reports analysis steps.
#'
#' @return An object of class \code{CEMiTool}.
#'
#' @examples
#' splitted_mods <- split_modules(cem_obj)
#' 
#' @rdname split_modules
#' @export
setGeneric('split_modules', function(cem_obj, ...) {
    standardGeneric('split_modules')
})

#' @rdname split_modules
setMethod('split_modules', signature(cem_obj='CEMiTool'),
          function(cem_obj, min_ngen=30, verbose=F) {

            if (verbose) {
                message('Splitting modules')
            }
            modules <- unique(cem_obj@module[, 'modules'])
            submods <- lapply(modules, function(mod){
        
                # subsets from exprs all genes inside module mod
                genes <- cem_obj@module[cem_obj@module[,'modules']==mod, 'genes']
                mod_exprs <- cem_obj@expression[genes,]
        
                # recalculates the correlation matrix for genes in module mod
                gene_cors <- cor(t(mod_exprs), use='everything', method='pearson')
        
                # checks if there are positive and negative correlations
                # inside module mod
                signs <- sign(range(gene_cors))
                if (signs[1] != signs[2]) {
                    # there are negative correlations inside module mod
        
                    if (verbose) {
                        message(paste0('Splitting module '), mod)
                    }
        
                    gene_dists <- as.dist(1 - gene_cors)
                    gene_clust <- hclust(gene_dists)
                    k <- cutree(gene_clust, 2)
                    pos_cors <- rownames(gene_cors)[which(k==1)]
                    neg_cors <- rownames(gene_cors)[which(k==2)]
        
                    # checks for the minimum module size of positive correlated genes
                    if (length(pos_cors) >= min_ngen) {
                        pos_mod <- data.frame(genes=pos_cors, modules=paste0(mod, '.A'))
                    } else {
                        pos_mod <- data.frame()
                    }
        
                    # checks for the minimum module size of negative correlated genes
                    if (length(neg_cors) >= min_ngen) {
                        neg_mod <- data.frame(genes=neg_cors, modules=paste0(mod, '.B'))
                    } else {
                        neg_mod <- data.frame()
                    }
        
                    splitted_mods <- rbind(pos_mod, neg_mod)
        
                    return(splitted_mods)
        
                } else {
                    # no negative correlations inside module mod
                    same_mod <- data.frame(genes=genes, modules=mod)
        
                    if (verbose) {
                        message(paste0('Module ', mod, ' will not be splitted'))
                    }
        
                    return(same_mod)
                }
            }) # end of lapply
            
            cem_obj@module <- do.call(rbind, submods)
            return(cem_obj)
        }
)


#' Co-expression module summarization 
#'
#' Summarizes modules using some statistics. 
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param method A character string indicating which summarization method 
#'                   is to be used. Default 'mean'. 
#' @param verbose Logical. If \code{TRUE}, reports analysis steps.
#'
#' @return A \code{data.frame} with summarized values.
#'
#'
#'
#' @examples
#' mod_summary <- mod_summary(cem_obj)
#'
#' @rdname mod_summary
#' @export
setGeneric('mod_summary', function(cem_obj, ...) {
    standardGeneric('mod_summary')
})

#' @rdname mod_summary
setMethod('mod_summary', signature(cem_obj='CEMiTool'),
          function(cem_obj, method=c('mean', 'eigengene'),
                   verbose=F)
          {
              method <- match.arg(method)

              if (verbose) {
                  message(paste0('Summarizing modules by ', method))
              }

              modules <- unique(cem_obj@module[, 'modules'])

              # mean expression of genes in modules
              if (method == 'mean') {
                  exprs <- data.table(cem_obj@expression, keep.rownames=T)
                  exprs_melt <- melt(exprs, id='rn', variable.name='samples',
                                     value.name='expression')
                  exprs_melt <- merge(exprs_melt, cem_obj@module, by.x='rn', by.y='genes')
                  summarized <- exprs_melt[, .(mean=mean(expression)),
                                           by=c('samples', 'modules')]
                  summarized <- dcast(summarized, modules~samples, value.var='mean')
                  setDF(summarized)

                  return(summarized)
                  # eigengene for each module
              } else if (method == 'eigengene') {
                  exprs_t <- t(cem_obj@expression)
                  colnames(exprs_t) <- rownames(exprs)
                  rownames(exprs_t) <- colnames(exprs)
                  me_list <- WGCNA::moduleEigengenes(exprs_t, 
                                                     colors=cem_obj@module[,2])
                  me_eigen <- data.table(t(me_list$eigengenes), keep.rownames=T)
                  setnames(me_eigen, c('modules', colnames(exprs)))
                  me_eigen[, modules := gsub('^ME', '', modules)]
                  setDF(me_eigen)

                  return(me_eigen)

              }
          }
)

