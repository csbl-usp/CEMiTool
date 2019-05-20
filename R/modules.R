#' @importFrom grDevices pdf
#' @importFrom grDevices colorRampPalette
#' @importFrom data.table data.table melt dcast setnames
#' @import WGCNA

.datatable.aware <- TRUE
#' Co-expression modules definition
#'
#' Defines co-expression modules
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param cor_method A character string indicating which correlation coefficient
#'                   is to be computed. Default \code{"pearson"}
#' @param cor_function A character string indicating the correlation function to be used. Default \code{'cor'}.
#' @param eps A value for accepted R-squared interval between subsequent beta values. Default is 0.1.
#' @param set_beta A value to override the automatically selected beta value. Default is NULL.
#' @param force_beta Whether or not to automatically force a beta value based on number of samples. Default is FALSE.
#' @param min_ngen Minimal number of genes per submodule. Default \code{20}.
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
#' data(expr0)
#' # Initialize CEMiTool object with expression
#' cem <- new_cem(expr0, filter=TRUE, apply_vst=FALSE)
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
                  cor_function='cor', eps=0.1,
                  set_beta=NULL, force_beta=FALSE,
                  min_ngen=20, merge_similar=TRUE,
                  network_type='unsigned', tom_type='signed',
                  diss_thresh=0.8, verbose=FALSE) {

    cor_method <- match.arg(cor_method)

    expr <- expr_data(cem, filter=cem@parameters$filter,
                      apply_vst=cem@parameters$apply_vst)
    if(nrow(expr) == 0){
        stop("CEMiTool object has no expression file!")
    }

    if(!is.null(set_beta) & force_beta){
        stop("Please specify only set_beta or force_beta!")
    }

    #vars <- mget(ls())
    #vars$expr <- NULL
    #cem <- get_args(cem, vars=vars)

    cem@parameters$cor_method <- cor_method
    cem@parameters$min_ngen <- min_ngen
    cem@parameters$merge_similar <- merge_similar
    cem@parameters$diss_thresh <- diss_thresh

    expr_t <- t(expr)
    names(expr_t) <- rownames(expr)
    rownames(expr_t) <- colnames(expr)

    if(cor_function == 'cor'){
        cor_options <- list(use="p", method=cor_method)
    }else if (cor_function == 'bicor'){
        cor_options <- list(use="p")
    }

    beta_data <- get_beta_data(cem, network_type, cor_function, cor_method, verbose)
    fit_indices <- beta_data$fitIndices
    wgcna_beta <- beta_data$powerEstimate

    if(!is.null(set_beta)){
        if(is.numeric(set_beta) & !(set_beta) %in% fit_indices[, "Power"]){
            stop(c("Parameter set_beta must be one of: ", paste(fit_indices$Power, collapse=" ")))
        }else if(is.character(set_beta) & tolower(set_beta)!= "wgcna"){
            stop("Unrecognized string character for parameter set_beta")
        }else if(is.na(set_beta)){
            stop("Parameter set_beta cannot be NA.")
        }
    }

    cem@fit_indices <- fit_indices

    k <- fit_indices[, 5]
    phi <- get_phi(cem)

    # Get CEMiTool Beta and R2
    if(is.null(set_beta)){
        r2_beta <- get_cemitool_r2_beta(cem, eps=eps)
        beta <- as.integer(r2_beta[2])
        r2 <- r2_beta[1]
    }else if(is.numeric(set_beta)){
        beta <- set_beta
        r2 <- fit_indices[fit_indices$Power == beta, 2]
    }else if(tolower(set_beta) == "wgcna"){
        beta <- wgcna_beta
        r2 <- fit_indices[fit_indices$Power == beta, 2]
    }

    if(force_beta){
        beta <- get_forced_beta(cem, network_type=network_type)
        r2 <- fit_indices[fit_indices$Power == beta, 2]
    }

    # Get network connectivities
    our_k <- get_connectivity(cem, beta)

    if(is.na(beta)){
        cem@parameters <- c(cem@parameters, NA)
        names(cem@parameters)[length(cem@parameters)] <- "beta"
        message('Could not specify the parameter Beta. No modules found.')
        return(cem)
    }

    cem@parameters$r2 <- r2
    cem@parameters$beta <- beta
    cem@parameters$phi <- phi

    # Get adjacency matrix
    cem <- get_adj(cem, beta=beta, network_type=network_type,
                   cor_function=cor_function, cor_method=cor_method)

    # Get modules
    mods <- get_mods(cem, tom_type=tom_type, min_ngen=min_ngen)
    n_mods <- length(unique(mods))
    cem@parameters$n_mods <- n_mods

    if(all(mods == 0)){
        message("No significant co-expression was found - all genes placed in the 'Not.Correlated' module")
        cem@parameters$n_mods <- NA
        return(cem)
    }

    # Number of modules
    if (n_mods == 0) {
        message('No modules found.')
        return(cem)
    }

    # if merge_similar=TRUE, merges similar modules
    if (merge_similar) {
        mods <- get_merged_mods(cem, mods, diss_thresh, verbose)
        # Re-count number of modules after merging
        n_mods <- length(unique(mods))
    }

    original_names <- setdiff(names(sort(table(mods), decreasing=TRUE)), 0)
    new_names <- paste0("M", 1:length(original_names))
    names(new_names) <- original_names
    new_names["0"] <- "Not.Correlated"

    out <- data.frame(genes=rownames(expr),
                      modules=new_names[as.character(mods)],
                      stringsAsFactors=FALSE)

    #params <- list('cor_method'=cor_method,
    #               'min_ngen'=min_ngen,
    #               'merge_similar'=merge_similar,
    #               'diss_thresh'=diss_thresh,
    #               'r2'=r2,
    #               'beta'=beta,
    #               'phi'=phi,
    #               'n_mods'=n_mods)

    #cem@parameters <- c(cem@parameters, params)
    cem@module <- out
    return(cem)

    })

#' Retrieve scale-free model fit data
#'
#' @param cem Object of class \code{CEMiTool}
#'
#' @return Object of class \code{data.frame}
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Get modules and beta data
#' cem <- find_modules(cem)
#' # Get fit data
#' fit_data(cem)
#'
#' @rdname fit_data
#' @export
setGeneric("fit_data", function(cem) {
    standardGeneric("fit_data")
})

#' @rdname fit_data
setMethod("fit_data", signature("CEMiTool"),
    function(cem){
        return(cem@fit_indices)
})

#' Soft-threshold beta data
#'
#' This function takes the input parameters from find_modules
#' and calculates the WGCNA soft-threshold parameters and returns
#' them.
#'
#' @param cem A CEMiTool object containing expression data
#' @param network_type A character string indicating to use either "unsigned"
#'           (default) or "signed" networks. Default \code{"unsigned"}.
#' @param cor_function A character string indicating the correlation function
#'           to be used. Default \code{'cor'}.
#' @param cor_method A character string indicating which correlation
#'           coefficient is to be computed. Default \code{"pearson"}
#' @param verbose Logical. If \code{TRUE}, reports analysis steps. Default \code{FALSE}
#' @param ... Optional parameters.
#'
#' @return A list containing the soft-threshold selected by WGCNA and scale-free model parameters
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr0, filter=TRUE, apply_vst=FALSE)
#' # Get beta data
#' beta_data <- get_beta_data(cem)
#'
#' @rdname get_beta_data
#' @export
setGeneric('get_beta_data', function(cem, ...) {
    standardGeneric('get_beta_data')
})
#' @rdname get_beta_data
#' @export
setMethod('get_beta_data', signature('CEMiTool'),
    function(cem, network_type="unsigned", cor_function="cor",
                      cor_method="pearson", verbose=FALSE){

    expr <- expr_data(cem, filter=cem@parameters$filter,
                      apply_vst=cem@parameters$apply_vst)
    if(nrow(expr) == 0){
        stop("CEMiTool object has no expression file!")
    }
    expr_t <- t(expr)
    names(expr_t) <- rownames(expr)
    rownames(expr_t) <- colnames(expr)

    if (verbose) {
        message('Selecting Beta')
        verbosity <- 10
    } else {
        verbosity <- 0
    }

    if(cor_function == 'cor'){
        cor_options <- list(use="p", method=cor_method)
    }else if (cor_function == 'bicor'){
        cor_options <- list(use="p")
    }

    # Define a range of soft-thresholding candidates
    if(network_type=="unsigned"){
        powers_end <- 20
    } else if (network_type=="signed"){
        powers_end <- 30
    }

    powers <- c(c(1:10), seq(12, powers_end, 2))

    ## Automatic selection of soft-thresholding power beta ##
    beta_data <- WGCNA::pickSoftThreshold(expr_t, powerVector=powers,
                                        networkType=network_type, moreNetworkConcepts=TRUE,
                                        corFnc=get('cor_function'),
                                        corOptions=cor_options,
                                        verbose=verbosity)
    return(beta_data)
})

#' Calculate phi
#'
#' This function takes a CEMiTool object and returns the phi parameter.
#'
#' @param cem A CEMiTool object containing the fit_indices slot
#' @param ... Optional parameters.
#'
#' @return The phi parameter
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr0, filter=TRUE, apply_vst=FALSE)
#' # Get modules and beta data
#' cem <- find_modules(cem)
#' # Get phi
#' get_phi(cem)
#'
#' @rdname get_phi
#' @export
setGeneric('get_phi', function(cem, ...) {
    standardGeneric('get_phi')
})
#' @rdname get_phi
#' @export
setMethod('get_phi', signature('CEMiTool'),
    function(cem){
        fit_indices <- fit_data(cem)
        if(nrow(fit_indices) == 0){
            stop("CEMiTool object has no fit_indices slot.")
        }
        fit <- -sign(fit_indices[, "slope"])*fit_indices[, "SFT.R.sq"]

        powers <- fit_indices$Power
        At  <- powers[length(powers)] - powers[1]
        A   <- 0.0
        A <- sum(0.5 * (tail(fit, -1) + head(fit, -1)) * (tail(powers, -1) - head(powers, -1)))

        # Area under the curve/threshold
        phi <- A/At
        return(phi)
})

#' Calculate CEMiTool beta and R2 values
#'
#' This function takes a CEMiTool object with beta data and returns
#' a vector containing the chosen beta and corresponding R squared value.
#'
#' @param cem A CEMiTool object containing the fit_indices slot
#' @param eps A value indicating the accepted interval between successive
#'                 values of R squared to use to calculate the selected beta.
#'              Default: 0.1.
#' @param ... Optional parameters.
#'
#' @return A vector containing R squared value and the chosen beta parameter.
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr0, filter=TRUE, apply_vst=FALSE)
#' # Get modules and beta data
#' cem <- find_modules(cem)
#' # Get CEMiTool R2 and beta values
#' get_cemitool_r2_beta(cem)
#'
#' @rdname get_cemitool_r2_beta
#' @export
setGeneric('get_cemitool_r2_beta', function(cem, ...) {
    standardGeneric('get_cemitool_r2_beta')
})
#' @rdname get_cemitool_r2_beta
setMethod('get_cemitool_r2_beta', signature(cem='CEMiTool'),
    function(cem, eps=0.1){
        fit_indices <- fit_data(cem)
        if(nrow(fit_indices) == 0){
            stop("CEMiTool object has no fit_indices slot.")
        }

        k <- fit_indices[, "mean.k."]
        powers <- fit_indices$Power
        fit <- -sign(fit_indices[, "slope"])*fit_indices[, "SFT.R.sq"]
        st <- c(NA,NA)

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
        return(st)
    })

#' Calculate network connectivity
#'
#' This function takes a CEMiTool object and returns the network connectivity.
#' @param cem Object of class \code{CEMiTool} containing the fit_indices slot
#' @param beta A soft-thresholding value to be used for the network.
#' @param ... Optional parameters.
#'
#' @return The value of the network's connectivity.
#'
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr0, filter=TRUE, apply_vst=FALSE)
#' # Get modules and beta data
#' cem <- find_modules(cem)
#' # Get network connectivity with example beta value 8
#' get_connectivity(cem, beta=8)
#'
#' @rdname get_connectivity
#' @export
setGeneric('get_connectivity', function(cem, ...) {
    standardGeneric('get_connectivity')
})
#' @rdname get_connectivity
setMethod('get_connectivity', signature(cem='CEMiTool'),
    function(cem, beta){
        fit_indices <- fit_data(cem)
        if(nrow(fit_indices) == 0){
            stop("CEMiTool object has no fit_indices slot.")
        }

        our_k <- NA
        if(!is.na(beta)) {
            if(beta <= 10) {
                our_k <- fit_indices[(beta) , 5]
            } else {
                line <- (beta + 10)/2
                our_k <- fit_indices[line, 5]
            }
        }
        return(our_k)
    })

#' Get or set adjacency matrix value
#'
#' This function takes a \code{CEMiTool} object containing expression values
#' and returns a CEMiTool object with an adjacency matrix in the adjacency slot.
#'
#' @param cem Object of class \code{CEMiTool}
#' @param value Object of class \code{matrix} containing adjacency data. Only used
#'           for setting adjacency values to CEMiTool object.
#' @param ... Optional parameters.
#'
#' @return Object of class \code{matrix} with adjacency values or object of class \code{CEMiTool}.
#'
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Initialize new CEMiTool object with expression
#' cem <- new_cem(expr0, filter=TRUE, apply_vst=FALSE)
#' # Calculate adjacency matrix with example beta value 8
#' cem <- get_adj(cem, beta=8)
#' # Return adjacency matrix
#' adj <- adj_data(cem)
#' # Check result
#' adj[1:5, 1:5]
#' # Set adjacency matrix to CEMiTool object
#' adj_data(cem) <- adj
#'
#' @rdname adj_data
#' @export
setGeneric('adj_data', function(cem, ...) {
    standardGeneric('adj_data')
})
#' @rdname adj_data
setMethod("adj_data", signature("CEMiTool"),
    function(cem) {
        return(cem@adjacency)
    })

#' @rdname adj_data
#' @export
setGeneric("adj_data<-", function(cem, value) {
    standardGeneric("adj_data<-")
})


#' @rdname adj_data
setReplaceMethod('adj_data', signature(cem='CEMiTool'),
    function(cem, value) {

        if(!is.matrix(value)){
            stop("The object provided is not a matrix object")
        }

        expr <- expr_data(cem, filter=cem@parameters$filter,
                          apply_vst=cem@parameters$apply_vst)
        if(nrow(expr) == 0){
            stop("CEMiTool object has no expression file!")
        }

        if(!identical(rownames(value), rownames(expr))){
            stop("Adjacency matrix provided does not reflect the names in the expression data.")
        }

        cem@adjacency <- value
         return(cem)
    })

#' Calculate adjacency matrix
#'
#' This function takes a \code{CEMiTool} object
#' and returns an adjacency matrix.
#'
#' @param cem Object of class \code{CEMiTool}
#' @param beta Selected soft-threshold value
#' @param network_type A character string indicating to use either "unsigned"
#'        (default) or "signed" networks. Default \code{"unsigned"}.
#' @param cor_function A character string indicating the correlation function
#'        to be used. Default \code{'cor'}.
#' @param cor_method A character string indicating which correlation
#'        coefficient is to be computed. Default \code{"pearson"}.
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with adjacency data
#'
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr0, filter=TRUE, apply_vst=FALSE)
#' # Calculate adjacency matrix with example beta value 8
#' cem <- get_adj(cem, beta=8)
#' # Check results
#' adj <- adj_data(cem)
#' adj[1:5, 1:5]
#'
#' @rdname get_adj
#' @export
setGeneric('get_adj', function(cem, ...) {
    standardGeneric('get_adj')
})
#' @rdname get_adj
setMethod("get_adj", signature("CEMiTool"),
    function(cem, beta, network_type="unsigned",
                cor_function="cor", cor_method="pearson") {

        if(missing(beta)){
            stop("Please provide a soft-threshold beta value. Run get_cemitool_r2_beta() for CEMiTool's default value.")
        }

        expr <- expr_data(cem, filter=cem@parameters$filter,
                          apply_vst=cem@parameters$apply_vst)
        if(nrow(expr) == 0){
            stop("CEMiTool object has no expression file!")
        }

        expr_t <- t(expr)
        names(expr_t) <- rownames(expr)
        rownames(expr_t) <- colnames(expr)

        if(cor_function == 'cor'){
            cor_options <- list(use="p", method=cor_method)
        }else if (cor_function == 'bicor'){
            cor_options <- list(use="p")
        }

        # Calculating adjacency matrix
        adj <- WGCNA::adjacency(expr_t, power=beta, type=network_type,
                                corFnc=cor_function, corOptions=cor_options)
        cem@adjacency <- adj
        return(cem)
    })

#' Calculate co-expression modules
#'
#' This function takes a \code{CEMiTool} object containing an adjacency matrix
#' together with the given network parameters, and returns the given
#' co-expression modules.
#' @param cem Object of class \code{CEMiTool}.
#' @param cor_function A character string indicating the correlation function
#'        to be used. Default \code{'cor'}.
#' @param cor_method A character string indicating which correlation
#'        coefficient is to be computed. Default \code{"pearson"}.
#' @param tom_type A character string indicating to use either "unsigned" or
#'           "signed" (default) TOM similarity measure.
#' @param min_ngen Minimal number of genes per module (Default: 20).
#' @param ... Optional parameters.
#'
#' @return Numeric labels assigning genes to modules.
#'
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr0, filter=TRUE, apply_vst=FALSE)
#' # Calculate adjacency matrix with example beta value 8
#' cem <- get_adj(cem, beta=8)
#' # Get module labels
#' mods <- get_mods(cem)

#' @rdname get_mods
#' @export
setGeneric('get_mods', function(cem, ...) {
    standardGeneric('get_mods')
})
#' @rdname get_mods
setMethod('get_mods', signature(cem='CEMiTool'),
    function(cem, cor_function="cor", cor_method="pearson",
             tom_type="signed", min_ngen=20) {

    expr <- expr_data(cem, filter=cem@parameters$filter,
                      apply_vst=cem@parameters$apply_vst)
    if(nrow(expr) == 0){
        stop("CEMiTool object has no expression file!")
    }

    adj <- adj_data(cem)
    if(nrow(adj) == 0){
        stop("CEMiTool object has no adjacency matrix!")
    }

    expr_t <- t(expr)
    names(expr_t) <- rownames(expr)
    rownames(expr_t) <- colnames(expr)

    if(cor_function == 'cor'){
        cor_options <- list(use="p", method=cor_method)
    }else if (cor_function == 'bicor'){
        cor_options <- list(use="p", method="pearson")
    }

    # Calculating Topological Overlap Matrix
    if (tom_type == 'signed') {
        tom <- WGCNA::TOMsimilarity(adj*sign(WGCNA::cor(expr_t, use=cor_options$use,
                                                        method=cor_options$method)), TOMType=tom_type)
    } else if (tom_type == 'unsigned') {
        tom <- WGCNA::TOMsimilarity(adj, TOMType=tom_type)
    }

    # Determining TOM based distance measure
    diss <- 1 - tom

    # Clustering
    tree <- hclust(as.dist(diss), method = 'average')

    # Cutting tree to determine modules
    mods <- dynamicTreeCut::cutreeDynamic(dendro = tree, distM = diss,
                                          deepSplit = 2,
                                          pamRespectsDendro = FALSE,
                                          minClusterSize = min_ngen)
    return(mods)
})

#' Merge similar modules
#'
#' This function takes a CEMiTool object with expression and a vector of
#' numeric labels to merge similar modules.
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param mods A vector containing numeric labels for each module gene
#' @param diss_thresh A threshold of dissimilarity for modules. Default is 0.8.
#' @param verbose Logical. If \code{TRUE}, reports analysis steps. Default \code{FALSE}
#' @param ... Optional parameters.
#'
#' @return Numeric labels assigning genes to modules.
#'
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr0, filter=TRUE, apply_vst=FALSE)
#' # Calculate adjacency matrix with example beta value 8
#' cem <- get_adj(cem, beta=8)
#' # Get modules
#' mods <- get_mods(cem)
#' # Merge similar modules
#' merged_mods <- get_merged_mods(cem, mods)
#'
#' @rdname get_merged_mods
#' @export
setGeneric('get_merged_mods', function(cem, ...) {
    standardGeneric('get_merged_mods')
})
#' @rdname get_merged_mods
setMethod('get_merged_mods', signature(cem='CEMiTool'),
    function(cem, mods, diss_thresh=0.8, verbose=FALSE){
        expr <- expr_data(cem, filter=cem@parameters$filter,
                          apply_vst=cem@parameters$apply_vst)
        if(nrow(expr) == 0){
            stop("CEMiTool object has no expression file!")
        }
        expr_t <- t(expr)
        names(expr_t) <- rownames(expr)
        rownames(expr_t) <- colnames(expr)

        if (verbose) {
            message('Merging modules based on eigengene similarity')
        }

        # Calculates eigengenes
        me_list <- WGCNA::moduleEigengenes(expr_t, colors=mods, grey=0)
        me_eigen <- me_list$eigengenes

        # Calculates dissimilarity of module eigengenes
        me_diss <- 1 - stats::cor(me_eigen)

        # Clustering module eigengenes
        me_tree <- hclust(as.dist(me_diss), method='average')

        # Setting cut height
        me_diss_thresh <- 1 - diss_thresh

        # Merging modules
        merged_mods <-  WGCNA::mergeCloseModules(expr_t, mods,
                                                 cutHeight=me_diss_thresh)

        # The merged modules colors
        mods <- merged_mods$colors

        return(mods)
    })

#' Co-expression module summarization
#'
#' Summarizes modules using mean or eigengene expression.
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param method A character string indicating which summarization method
#' is to be used. Can be 'eigengene', 'mean' or 'median'. Default is 'mean'.
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
#' @rdname mod_summary
#' @export
setGeneric('mod_summary', function(cem, ...) {
    standardGeneric('mod_summary')
})

#' @rdname mod_summary
setMethod('mod_summary', signature(cem='CEMiTool'),
    function(cem, method=c('mean', 'median', 'eigengene'),
             verbose=FALSE){

        method <- match.arg(method)

        if(length(cem@module) == 0){
            stop("No modules in CEMiTool object! Did you run find_modules()?")
          }

        if (verbose) {
            message(paste0('Summarizing modules by ', method))
        }

        modules <- unique(cem@module[, 'modules'])

        expr <- expr_data(cem, filter=cem@parameters$filter,
                          apply_vst=cem@parameters$apply_vst)
        if(nrow(expr) == 0){
            stop("CEMiTool object has no expression file!")
        }

        # mean expression of genes in modules
        if (method == 'mean' | method == 'median') {
            func <- get(method)
            expr <- data.table(expr, keep.rownames=TRUE)
            expr_melt <- data.table::melt(expr, id='rn', variable.name='samples',
                               value.name='expression')
            expr_melt <- merge(expr_melt, cem@module, by.x='rn', by.y='genes')
            summarized <- expr_melt[, list(method=as.double(func(expression))),
                                     by=c('samples', 'modules')]
            summarized <- data.table::dcast(summarized, modules~samples, value.var="method")
            data.table::setDF(summarized)

            return(summarized)
            # eigengene for each module
        } else if (method == 'eigengene') {
            expr_t <- t(expr)
            colnames(expr_t) <- rownames(expr)
            rownames(expr_t) <- colnames(expr)
            me_list <- WGCNA::moduleEigengenes(expr_t,
                                               colors=cem@module[,2])
            me_eigen <- data.table(t(me_list$eigengenes), keep.rownames=TRUE)
            setnames(me_eigen, c('modules', colnames(expr)))
            me_eigen[, modules := gsub('^ME', '', modules)]
            data.table::setDF(me_eigen)

            return(me_eigen)
        }
    })

#' Get hubs
#'
#' Returns \code{n} genes in each module with high connectivity.
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param n Number of hubs to return from each module. If "all", returns all genes
#' in decreasing order of connectivity. Default: 5.
#' @param method Method for hub calculation. Either "adjacency" or "kME".
#' Default: "adjacency"
#' @param ... Optional parameters.
#'
#' @return A \code{list} containing hub genes for each module and the value of
#' the calculated method.
#'
#' @examples
#' # Get example CEMiTool object
#' data(cem)
#' # Get module hubs
#' hubs <- get_hubs(cem, n=10, "adjacency")
#'
#' @rdname get_hubs
#' @export
setGeneric('get_hubs', function(cem, ...) {
    standardGeneric('get_hubs')
})

#' @rdname get_hubs
setMethod('get_hubs', signature(cem='CEMiTool'),
      function(cem, n=5, method="adjacency"){
          if(nrow(cem@adjacency) == 0){
              stop("Make sure that you ran the method find_modules.")
          }

          if(method == "adjacency"){
              mod2gene <- split(cem@module$genes, cem@module$modules)
              hubs <- lapply(mod2gene, function(x){
                  if (length(x) > 1) {
                      mod_adj <- cem@adjacency[x, x]
                      diag(mod_adj) <- NA
                      if(n == "all"){
                          top <- sort(rowSums(mod_adj, na.rm=TRUE), decreasing=TRUE)
                      }else{
                          top <- head(sort(rowSums(mod_adj, na.rm=TRUE), decreasing=TRUE), n=n)
                      }
                      return(top)
                  } else {
                      return(character())
                  }
              })
              return(hubs)
          }else if(method == "kME"){
              eigens <- mod_summary(cem, "eigengene")
              rownames(eigens) <- eigens$modules
              eigens$modules <- NULL
              eigens <- as.data.frame(t(eigens))

              expr <- expr_data(cem, filter=cem@parameters$filter,
                                apply_vst=cem@parameters$apply_vst)
              datExpr <- as.data.frame(t(expr))
              kmes <- signedKME(datExpr, eigens)
              names(kmes) <- names(eigens)

              kme_list <- as.list(kmes)
              kme_list <- lapply(kme_list, function(x){
                  names(x) <- rownames(kmes)
                  x
              })
              hubs <- lapply(kme_list, function(x){
                  if(n == "all"){
                      top <- sort(x, decreasing=TRUE)
                  }else{
                      top <- head(sort(x, decreasing=TRUE), n = n)
                  }
                  return(top)
              })
              return(hubs)
          }
    })
