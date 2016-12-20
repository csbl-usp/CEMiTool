library(WGCNA)

#' Co-expression modules definition
#'
#' Defines co-expression modules
#'
#' @param exprs gene expression \code{data.frame}
#' @param cor_method a character string indicating which correlation coefficient is to be computed
#' @param verbose logical. Report analysis steps
#'
#' @return just god knows 
#'
#' @examples
#' coex_mods <- find_modules(exprs=expression.df, cor_method='pearson')
#'
#' @export
find_modules <- function(exprs, cor_method=c('pearson', 'spearman'),
                             min_mod_size=30, merge_similar=T,
                             diss_thresh=0.8, verbose=F)
{
    
    if (is.null(exprs)){
        stop('Must provide expression data')
    }

    exprs_t <- t(exprs)
    names(exprs_t) <- rownames(exprs)
    rownames(exprs_t) <- colnames(exprs)


    if (verbose) {
        message('Selecting Beta')
    }

    # Define a range of soft-thresholding candidates
    powers <- c(c(1:10), seq(12, 20, 2))
    
    ## Automatic selection of soft-thresholding power beta ##
    beta <- pickSoftThreshold(exprs_t, powerVector=powers, verbose=5,
                              networkType='signed', moreNetworkConcepts=TRUE,
                              corOptions=list(use='p',
                                              method=match.arg(cor_method)))
    
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
    our_adj <- adjacency(exprs_t, power=our_beta, type='signed')

    # Calculating Topological Overlap Matrix
    our_tom <- TOMsimilarity(our_adj)

    # Determining TOM based distance measure
    our_diss <- 1 - our_tom

    # Clustering
    our_tree <- hclust(as.dist(our_diss), method = 'average')

    # Cutting tree to determine modules
    our_mods <- cutreeDynamic(dendro = our_tree, distM = our_diss,
                              deepSplit = 2,
                              pamRespectsDendro = FALSE,
                              minClusterSize = min_mod_size)
    our_table <- table(our_mods)

    our_colors <- labels2colors(our_mods)

    out <- data.frame(genes=rownames(exprs), modules=our_colors)

    n_mods <- length(unique(our_colors))

    if (n_mods <= 1) {
        stop('Could not specify the parameter Beta. No modules found.')
    }

    # if merge_similar=TRUE, merges similar modules
    if (merge_similar) {
        if (verbose) {
            message('Merging modules based on eigengene similarity')
        }
        # Calcultes eigengenes
        me_list <- moduleEigengenes(exprs_t, colors=our_colors)
        me_eigen <- me_list$eigengenes

        # Calculates dissimilarity of module eigengenes
        me_diss <- 1 - cor(me_eigen)

        # Clustering module eigengenes
        me_tree <- hclust(as.dist(me_diss), method='average')

        # Setting cut height
        me_diss_thresh <- 1 - diss_thresh 

        # Merging modules
        merged_mods <-  mergeCloseModules(exprs_t, our_colors,
                                          cutHeight=diss_thresh)

        # The merged modules colors
        our_colors <- merged_mods$colors

        # Eigengenes of the new merged modules
        merged_mods <- merged_mods$newMEs

        out[, 'modules'] <- our_colors
    }

    return(out)
}



#' Co-expression module refinement 
#'
#' Refines modules by splitting them into submodules based on correlation sign.
#'
#' @param exprs gene expression \code{data.frame}
#' @param gene_module two column \code{data.frame}. First column with
#'        gene identifiers and second column with module information
#' @param verbose logical. Report analysis steps
#'
#' @return A \code{data.frame} with gene identifier and module information.
#'
#' @examples
#' splitted_mods <- split_modules(exprs=expression.df, gene_module=coex)
#'
#' @export
split_modules <- function(exprs, gene_module, verbose=F) {}



#' Co-expression module refinement 
#'
#' Refines modules by excluding non-significant edges. 
#'
#' @param exprs gene expression \code{data.frame}
#' @param gene_module two column \code{data.frame}. First column with
#'        gene identifiers and second column with module information
#' @param nperm numeric. Number of permutations for edge significance test
#' @param use_lpc logical. If TRUE uses Local Partial Correlation method for
#'        module refinement
#' @param verbose logical. Report analysis steps
#'
#' @return A list of with the following components:
#'
#'gene_module: \code{data.frame} with gene identifier and module information
#'edge_list: \code{data.frame} with the genes which consist of each edge in 
#'            a module.
#'
#'
#' @examples
#' refined_mods <- refine_modules(exprs=expression.df, gene_module=coex)
#'
#' @export
refine_modules <- function(exprs, gene_module, nperm=1000,
                         use_lpc=F, verbose=F){}


