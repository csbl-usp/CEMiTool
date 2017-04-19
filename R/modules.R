#' @importFrom grDevices pdf
#' @importFrom grDevices colorRampPalette
#' @import stats
#' @import data.table

.datatable.aware=TRUE
#' Co-expression modules definition
#'
#' Defines co-expression modules
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param cor_method A character string indicating which correlation coefficient
#'                   is to be computed.
#' @param min_ngen Minimal number of genes per submodule. Default \code{30}.
#' @param merge_similar Logical. If \code{TRUE}, merge similar modules.
#' @param network_type A character string indicating to use either "unsigned" (default) or "signed" networks.
#' @param tom_type A character string indicating to use either "unsigned" or "signed" (default) TOM similarity measure.
#' @param diss_thresh Module merging correlation threshold for eigengene similarity.
#'        Default \code{0.8}.
#' @param verbose Logical. If \code{TRUE}, reports analysis steps.
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} 
#'
#' @examples
#' cem <- new('CEMiTool', expression=expr)
#' cem <- find_modules(cem)
#'
#' @rdname find_modules
#' @export
setGeneric('find_modules', function(cem, ...) {
    standardGeneric('find_modules')
})

#' @rdname find_modules
#' @export
setMethod('find_modules', signature('CEMiTool'), 
          function(cem, cor_method=c('pearson', 'spearman'),
                             min_ngen=20, merge_similar=T,
                             network_type='unsigned', tom_type='signed',
                             diss_thresh=0.8, verbose=F)
{
    cor_method <- match.arg(cor_method)
    
    expr <- expr_data(cem)

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
    if(network_type=="unsigned"){
        powers_end <- 20 
    } else if (network_type=="signed"){
        powers_end <- 30
    }
    
    powers <- c(c(1:10), seq(12, powers_end, 2))
    
    ## Automatic selection of soft-thresholding power beta ##
    beta <- WGCNA::pickSoftThreshold(expr_t, powerVector=powers,
                              networkType=network_type, moreNetworkConcepts=TRUE,
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
    if(is.na(our_beta)){
        cem@parameters <- c(cem@parameters, NA)
        stop('Could not specify the parameter Beta. No modules found.')
    }
    
    our_r2 <- st[1]

    # Calculating adjacency matrix
    our_adj <- WGCNA::adjacency(expr_t, power=our_beta, type=network_type)
    cem@adjacency <- our_adj

    # Calculating Topological Overlap Matrix
    if (tom_type == 'signed') {
        our_tom <- WGCNA::TOMsimilarity(our_adj*sign(cor(expr_t)), TOMType=tom_type)
    } else if (tom_type == 'unsigned') {
        our_tom <- WGCNA::TOMsimilarity(our_adj, TOMType=tom_type)
    }

    # Determining TOM based distance measure
    our_diss <- 1 - our_tom

    # Clustering
    our_tree <- hclust(as.dist(our_diss), method = 'average')

    # Cutting tree to determine modules
    our_mods <- dynamicTreeCut::cutreeDynamic(dendro = our_tree, distM = our_diss,
                              deepSplit = 2,
                              pamRespectsDendro = FALSE,
                              minClusterSize = min_ngen)
    
    # Number of modules
    n_mods <- length(unique(our_mods))
    

    # checks number of modules found
    if (n_mods <= 1) {
        stop('No modules found.')
    }

    # if merge_similar=TRUE, merges similar modules
    if (merge_similar) {
        if (verbose) {
            message('Merging modules based on eigengene similarity')
        }
        # Calculates eigengenes
        me_list <- WGCNA::moduleEigengenes(expr_t, colors=our_mods, grey=0)
        me_eigen <- me_list$eigengenes

        # Calculates dissimilarity of module eigengenes
        me_diss <- 1 - cor(me_eigen)

        # Clustering module eigengenes
        me_tree <- hclust(as.dist(me_diss), method='average')

        # Setting cut height
        me_diss_thresh <- 1 - diss_thresh 

        # Merging modules

        merged_mods <-  WGCNA::mergeCloseModules(expr_t, our_mods,
                                          cutHeight=me_diss_thresh)

        # The merged modules colors
        our_mods <- merged_mods$colors
        
        # Number of modules after merging
        n_mods <- length(unique(our_mods))
    }

    original_names <- setdiff(names(sort(table(our_mods), decreasing=TRUE)), 0)
    new_names <- paste0("M", 1:length(original_names))
    names(new_names) <- original_names
    new_names["0"] <- "Not.Correlated"

    out <- data.frame(genes=rownames(expr), 
                      modules=new_names[as.character(our_mods)], 
                      stringsAsFactors=FALSE)
    
    params <- list('cor_method'=cor_method,
                   'min_ngen'=min_ngen,
                   'merge_similar'=merge_similar,
                   'diss_thresh'=diss_thresh,
                   'r2'=our_r2, 
                   'beta'=our_beta, 
                   'phi'=phi,
                   'n_mods'=n_mods
                   )
    cem@parameters <- c(cem@parameters, params)
    cem@module <- out
    cem@fit_indices <- beta$fitIndices
    return(cem)
})


#' Co-expression module refinement 
#'
#' Refines modules by splitting them into submodules based on correlation sign.
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param min_ngen Minimal number of genes per submodule. Default \code{30}.
#' @param verbose Logical. If \code{TRUE}, reports analysis steps.
#' @param ... Optional parameters.
#'
#' @return An object of class \code{CEMiTool}.
#'
#' @examples
#' splitted_mods <- split_modules(cem)
#' 
#' @rdname split_modules
#' @export
setGeneric('split_modules', function(cem, ...) {
    standardGeneric('split_modules')
})

#' @rdname split_modules
setMethod('split_modules', signature(cem='CEMiTool'),
          function(cem, min_ngen=30, verbose=F) {
            if (verbose) {
                message('Splitting modules')
            }
            modules <- unique(cem@module[, 'modules'])
            submods <- lapply(modules, function(mod){
                if (verbose) {
                    message('Evaluating Module ', mod)
                }
        
                # subsets from expr all genes inside module mod
                genes <- cem@module[cem@module[,'modules']==mod, 'genes']
                mod_expr <- expr_data(cem)[genes,]
        
                # recalculates the correlation matrix for genes in module mod
                gene_cors <- cor(t(mod_expr), use='everything', method='pearson')

                test_obj <- diptest::dip.test(gene_cors[lower.tri(gene_cors)])
                bimod_pval <- test_obj$p.value

                if(verbose){
                    tryCatch({
                        bimod_ratio <- modes::bimodality_ratio(gene_cors[lower.tri(gene_cors)])
                        message("Bimodality ratio: ", bimod_ratio)
                        message("Hartingans Dip Test - p-value: ", bimod_pval)
                    })
                }
                
        
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
        
                    # checks for the minimum module size of positive and negative correlated genes
                    if (length(pos_cors) >= min_ngen && length(neg_cors) >= min_ngen) {
                        modA_name <- paste0(mod, ".A")
                        modB_name <- paste0(mod, ".B")
                        if(verbose){
                            prop <- sum(gene_cors < 0, na.rm=T)/(nrow(gene_cors)**2 - nrow(gene_cors))
                            prop_pos <- sum(gene_cors[pos_cors, pos_cors] < 0, na.rm=T)/(length(pos_cors)**2 - length(pos_cors))
                            prop_neg <- sum(gene_cors[neg_cors, neg_cors] < 0, na.rm=T)/(length(neg_cors)**2 - length(neg_cors))
                            message("Proportion of negative correlations on ", mod, ": ", round(prop, digits=4))
                            message("Proportion of negative correlations on ", modA_name, ": ", round(prop_pos, digits=4))
                            message("Proportion of negative correlations on ", modB_name, ": ", round(prop_neg, digits=4))
                        }
                        pos_mod <- data.frame(genes=pos_cors, modules=modA_name, stringsAsFactors=FALSE)
                        neg_mod <- data.frame(genes=neg_cors, modules=modB_name, stringsAsFactors=FALSE)
                        mod_return <- rbind(pos_mod, neg_mod)
                    } else {
                        if (verbose) {
                            message(paste0('Module ', mod, ' will not be splitted'))
                        }
                        mod_return <- data.frame(genes=genes, modules=mod, stringsAsFactors=FALSE)
                    } 
                    return(mod_return)
        
                } else {
                    # no negative correlations inside module mod
                    same_mod <- data.frame(genes=genes, modules=mod, stringsAsFactors=FALSE)
        
                    if (verbose) {
                        message(paste0('Module ', mod, ' will not be splitted'))
                    }
        
                    return(same_mod)
                }
            }) # end of lapply
            
            cem@module <- do.call(rbind, submods)
            return(cem)
        }
)


#' Co-expression module summarization 
#'
#' Summarizes modules using some statistics. 
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param method A character string indicating which summarization method 
#'                   is to be used. Default 'mean'. 
#' @param verbose Logical. If \code{TRUE}, reports analysis steps.
#' @param ... Optional parameters.
#'
#' @return A \code{data.frame} with summarized values.
#'
#'
#'
#' @examples
#' mod_summary <- mod_summary(cem)
#'
#' @rdname mod_summary
#' @export
setGeneric('mod_summary', function(cem, ...) {
    standardGeneric('mod_summary')
})

#' @rdname mod_summary
setMethod('mod_summary', signature(cem='CEMiTool'),
          function(cem, method=c('mean', 'eigengene'),
                   verbose=F)
          {
              method <- match.arg(method)

              if (verbose) {
                  message(paste0('Summarizing modules by ', method))
              }

              modules <- unique(cem@module[, 'modules'])

              # mean expression of genes in modules
              if (method == 'mean') {
                  expr <- data.table(expr_data(cem), keep.rownames=T)
                  expr_melt <- melt(expr, id='rn', variable.name='samples',
                                     value.name='expression')
                  expr_melt <- merge(expr_melt, cem@module, by.x='rn', by.y='genes')
                  summarized <- expr_melt[, .(mean=mean(expression)),
                                           by=c('samples', 'modules')]
                  summarized <- dcast(summarized, modules~samples, value.var='mean')
                  setDF(summarized)

                  return(summarized)
                  # eigengene for each module
              } else if (method == 'eigengene') {
                  expr_t <- t(expr_data(cem))
                  colnames(expr_t) <- rownames(expr)
                  rownames(expr_t) <- colnames(expr)
                  me_list <- WGCNA::moduleEigengenes(expr_t, 
                                                     colors=cem@module[,2])
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
#' @param cem Object of class \code{CEMiTool}.
#' @param n Number of genes to return in each module (default: 5).
#' @param ... Optional parameters.
#'
#' @return A \code{list} containing hub genes.
#'
#'
#'
#' @examples
#'
#' @rdname get_hubs
#' @export
setGeneric('get_hubs', function(cem, ...) {
    standardGeneric('get_hubs')
})

#' @rdname get_hubs
setMethod('get_hubs', signature(cem='CEMiTool'),
          function(cem, n=5){
              if(nrow(cem@adjacency) == 0){
                  stop("Make sure that you ran the method find_modules.")
              }
              mod2gene <- split(cem@module$genes, cem@module$modules)
              hubs <- lapply(mod2gene, function(x){
                                           if (length(x) > 1) {
                                               mod_adj <- cem@adjacency[x, x]
                                               top <- head(sort(rowSums(mod_adj), decreasing=T), n=n)
                                               return(names(top))
                                           } else {
                                               return(character())
                                           }
                                       })
              return(hubs)
          })
