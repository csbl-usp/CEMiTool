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
#' @param ... other optional parameters
#'
#' @return Object of class \code{CEMiTool} 
#'
#' @examples
#' data(expr)
#' cem_obj <- new('CEMiTool', expression=expr)
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
    cor_method <- match.arg(cor_method)
    
    expr <- expr_data(cem_obj)

    expr_t <- t(expr)
    names(expr_t) <- rownames(expr)
    rownames(expr_t) <- colnames(expr)


    if (verbose) {
        message('Selecting Beta')
        verbosity <- 10
    } else {
        verbosity <- 0
    }

    # Define a range of soft-thresholding candidates
    powers <- c(c(1:10), seq(12, 20, 2))
    
    ## Automatic selection of soft-thresholding power beta ##
    beta <- WGCNA::pickSoftThreshold(expr_t, powerVector=powers,
                              networkType='signed', moreNetworkConcepts=TRUE,
                              corOptions=list(use='p',
                                              method=cor_method),
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
                j <- which.max(k[count:(count+2)]) + count - 1
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
    our_adj <- WGCNA::adjacency(expr_t, power=our_beta, type='signed')
	cem_obj@adjacency <- our_adj

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
    
    out <- data.frame(genes=rownames(expr), modules=our_colors, stringsAsFactors=FALSE)

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
        me_list <- WGCNA::moduleEigengenes(expr_t, colors=our_colors)
        me_eigen <- me_list$eigengenes

        # Calculates dissimilarity of module eigengenes
        me_diss <- 1 - cor(me_eigen)

        # Clustering module eigengenes
        me_tree <- hclust(as.dist(me_diss), method='average')

        # Setting cut height
        me_diss_thresh <- 1 - diss_thresh 

        # Merging modules
        merged_mods <-  WGCNA::mergeCloseModules(expr_t, our_colors,
                                          cutHeight=diss_thresh)

        # The merged modules colors
        merged_mods <- factor(merged_mods$colors)

        levels(merged_mods) <- paste0('M', rank(-table(merged_mods)))

        out[, 'modules'] <- as.character(merged_mods)
    }
    
    params <- list('cor_method'=cor_method,
                   'min_ngen'=min_ngen,
                   'merge_similar'=merge_similar,
                   'diss_thresh'=diss_thresh
                   )
    cem_obj@parameters <- append(params, list('r2'=our_r2, 'beta'=our_beta, 'phi'=phi))
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
#' @param ... other optional parameters
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
        
                # subsets from expr all genes inside module mod
                genes <- cem_obj@module[cem_obj@module[,'modules']==mod, 'genes']
                mod_expr <- expr_data(cem_obj)[genes,]
        
                # recalculates the correlation matrix for genes in module mod
                gene_cors <- cor(t(mod_expr), use='everything', method='pearson')
        
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
                        pos_mod <- data.frame(genes=pos_cors, modules=paste0(mod, '.A'), stringsAsFactors=FALSE)
                    } else {
                        pos_mod <- data.frame(stringsAsFactors=FALSE)

                    }
        
                    # checks for the minimum module size of negative correlated genes
                    if (length(neg_cors) >= min_ngen ) {
                        neg_mod <- data.frame(genes=neg_cors, modules=paste0(mod, '.B'), stringsAsFactors=FALSE)
                    } else {
                        neg_mod <- data.frame(stringsAsFactors=FALSE)
                    }
        
                    splitted_mods <- rbind(pos_mod, neg_mod)
        
                    return(splitted_mods)
        
                } else {
                    # no negative correlations inside module mod
                    same_mod <- data.frame(genes=genes, modules=mod, stringsAsFactors=FALSE)
        
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
#' @param ... other optional parameters
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
                  expr <- data.table(expr_data(cem_obj), keep.rownames=T)
                  expr_melt <- melt(expr, id='rn', variable.name='samples',
                                     value.name='expression')
                  expr_melt <- merge(expr_melt, cem_obj@module, by.x='rn', by.y='genes')
                  summarized <- expr_melt[, .(mean=mean(expression)),
                                           by=c('samples', 'modules')]
                  summarized <- dcast(summarized, modules~samples, value.var='mean')
                  setDF(summarized)

                  return(summarized)
                  # eigengene for each module
              } else if (method == 'eigengene') {
                  expr_t <- t(expr_data(cem_obj))
                  colnames(expr_t) <- rownames(expr)
                  rownames(expr_t) <- colnames(expr)
                  me_list <- WGCNA::moduleEigengenes(expr_t, 
                                                     colors=cem_obj@module[,2])
                  me_eigen <- data.table(t(me_list$eigengenes), keep.rownames=T)
                  setnames(me_eigen, c('modules', colnames(expr)))
                  me_eigen[, modules := gsub('^ME', '', modules)]
                  setDF(me_eigen)

                  return(me_eigen)

              }
          }
)


#' Get Hubs 
#'
#' Returns \code{n} genes in each module with high connectivity. 
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param n Number of genes to return in each module (default: 5).
#'
#' @return A \code{list} containing hub genes.
#'
#'
#'
#' @examples
#'
#' @rdname get_hubs
#' @export
setGeneric('get_hubs', function(cem_obj, ...) {
    standardGeneric('get_hubs')
})

#' @rdname get_hubs
setMethod('get_hubs', signature(cem_obj='CEMiTool'),
          function(cem_obj, n=5){
              if(nrow(cem_obj@adjacency) == 0){
                  stop("Make sure that you ran the method find_modules.")
              }
              mod2gene <- split(cem_obj@module$genes, cem_obj@module$modules)
              hubs <- lapply(mod2gene, function(x){
                                           mod_adj <- cem_obj@adjacency[x, x]
                                           top <- head(sort(rowSums(mod_adj), decreasing=T), n=n)
                                           return(names(top))
                                       })
              return(hubs)
          })
