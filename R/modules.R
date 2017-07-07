#' @importFrom grDevices pdf
#' @importFrom grDevices colorRampPalette
#' @importFrom data.table fread setDF
#' @import WGCNA

.datatable.aware=TRUE
#' Co-expression modules definition
#'
#' Defines co-expression modules
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param cor_method A character string indicating which correlation coefficient
#'                   is to be computed. Default \code{"pearson"}
#' @param cor_function A character string indicating the correlation function to be used. Default \code{'cor'}
#' @param min_ngen Minimal number of genes per submodule. Default \code{30}.
#' @param merge_similar Logical. If \code{TRUE}, (the default) merge similar modules.
#' @param network_type A character string indicating to use either "unsigned" (default) or "signed" networks. Default \code{"unsigned"}
#' @param tom_type A character string indicating to use either "unsigned" or "signed" (default) TOM similarity measure.
#' @param diss_thresh Module merging correlation threshold for eigengene similarity.
#'        Default \code{0.8}.
#' @param verbose Logical. If \code{TRUE}, reports analysis steps. Default \code{FALSE}
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} 
#'
#' @examples
#' # Get example expression data
#' data(expr)
#' # Initialize CEMiTool object with expression
#' cem <- new('CEMiTool', expression=expr)
#' # Filter data
#' cem <- filter_expr(cem)
#' # Define network modules
#' cem <- find_modules(cem)
#' # Check results
#' head(module_genes(cem))
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
		   	     cor_function='cor',
                             min_ngen=20, merge_similar=TRUE,
                             network_type='unsigned', tom_type='signed',
                             diss_thresh=0.8, verbose=FALSE)
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

    if(cor_function == 'cor'){
	cor_options <- list(use="p", method=cor_method)
    }else if (cor_function == 'bicor'){
	cor_options <- list(use="p")
    }
    
    ## Automatic selection of soft-thresholding power beta ##
    beta <- WGCNA::pickSoftThreshold(expr_t, powerVector=powers,
                              networkType=network_type, moreNetworkConcepts=TRUE,
			      corFnc=get('cor_function'),
                              corOptions=cor_options,
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
        our_tom <- WGCNA::TOMsimilarity(our_adj*sign(WGCNA::cor(expr_t)), TOMType=tom_type)
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
        me_diss <- 1 - stats::cor(me_eigen)

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

#' Co-expression module summarization 
#'
#' Summarizes modules using mean or eigengene expression.
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param method A character string indicating which summarization method 
#'                   is to be used. Default 'mean'. 
#' @param verbose Logical. If \code{TRUE}, reports analysis steps.
#' @param ... Optional parameters.
#'
#' @return A \code{data.frame} with summarized values.
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Summarize results
#' mod_summary <- mod_summary(cem)
#' 
#' 
#' @rdname mod_summary
#' @export
setGeneric('mod_summary', function(cem, ...) {
    standardGeneric('mod_summary')
})

#' @rdname mod_summary
setMethod('mod_summary', signature(cem='CEMiTool'),
          function(cem, method=c('mean', 'eigengene'),
                   verbose=FALSE)
          {
              method <- match.arg(method)

              if (verbose) {
                  message(paste0('Summarizing modules by ', method))
              }

              modules <- unique(cem@module[, 'modules'])

              # mean expression of genes in modules
              if (method == 'mean') {
                  expr <- data.table(expr_data(cem), keep.rownames=TRUE)
                  expr_melt <- melt(expr, id='rn', variable.name='samples',
                                     value.name='expression')
                  expr_melt <- merge(expr_melt, cem@module, by.x='rn', by.y='genes')
                  summarized <- expr_melt[, list(mean=mean(expression)),
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
                  me_eigen <- data.table(t(me_list$eigengenes), keep.rownames=TRUE)
                  setnames(me_eigen, c('modules', colnames(expr)))
                  me_eigen[, modules := gsub('^ME', '', modules)]
                  setDF(me_eigen)

                  return(me_eigen)

              }
          }
)


#' Get hubs 
#'
#' Returns \code{n} genes in each module with high connectivity. 
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param n Number of genes to return in each module (default: 5).
#' @param ... Optional parameters.
#'
#' @return A \code{list} containing hub genes.
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Get module hubs
#' hubs <- get_hubs(cem, n=10)
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
                                               top <- head(sort(rowSums(mod_adj), decreasing=TRUE), n=n)
                                               return(names(top))
                                           } else {
                                               return(character())
                                           }
                                       })
              return(hubs)
          })
