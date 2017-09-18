#' @importFrom grDevices pdf
#' @importFrom grDevices colorRampPalette
#' @import data.table
#' @import WGCNA

.datatable.aware <- TRUE
#' Co-expression modules definition
#'
#' Defines co-expression modules
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param cor_method A character string indicating which correlation coefficient
#'                   is to be computed. Default \code{"pearson"}
#' @param cor_function A character string indicating the correlation function to be used. Default \code{'cor'}
#' @param set_beta A value to override the automatically selected beta value. Default is NULL.
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
#' data(expr)
#' # Initialize CEMiTool object with expression
#' cem <- new_cem(expr)
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
                          cor_function='cor', set_beta=NULL,
                          min_ngen=20, merge_similar=TRUE,
                          network_type='unsigned', tom_type='signed',
                          diss_thresh=0.8, verbose=FALSE) {

	cor_method <- match.arg(cor_method)
    
    expr <- expr_data(cem)
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

	beta_data <- get_beta_data(cem, network_type, cor_function, cor_method, verbose)
    fit_indices <- beta_data$fitIndices
    wgcna_beta <- beta_data$powerEstimate
		    
	if(!is.null(set_beta)){
	    if(is.numeric(set_beta) & !(set_beta) %in% fit_indices[, 2]){
		    stop(c("Parameter set_beta must be one of: ", paste(fit_indices$Power, collapse=" ")))
	    }else if(is.character(set_beta) & tolower(set_beta)!= "wgcna"){
			stop("Unrecognized string character for parameter set_beta")
		}
	}	    
		    
	cem@fit_indices <- fit_indices
			    
	k <- fit_indices[, 5]
	phi <- get_phi(cem)

	# Get CEMiTool Beta and R2
    if(is.null(set_beta)){
	    r2_beta <- get_cemitool_r2_beta(cem, eps=0.1)    
		beta <- as.integer(r2_beta[2])
		r2 <- r2_beta[1]
	}else if(is.numeric(set_beta)){
		beta <- set_beta
		r2 <- fit_indices[fit_indices$Power == beta, 2]
	}else if(tolower(set_beta) == "wgcna"){
		beta <- wgcna_beta
		r2 <- fit_indices[fit_indices$Power == beta, 2]
	}
	    
	# Get network connectivities
	our_k <- get_connectivity(cem, beta)
		    
	if(is.na(beta)){
		cem@parameters <- c(cem@parameters, NA)
		stop('Could not specify the parameter Beta. No modules found.')
	}

	# Get modules
	mods <- get_mods(cem, beta=beta, network_type=network_type, 
					 tom_type=tom_type, cor_function=cor_function, 
					 cor_method=cor_method, min_ngen=min_ngen)
	    
	# Number of modules
	n_mods <- length(unique(mods))
	if (n_mods <= 1) {
		stop('No modules found.')
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
			    
	params <- list('cor_method'=cor_method,
	               'min_ngen'=min_ngen,
	               'merge_similar'=merge_similar,
	               'diss_thresh'=diss_thresh,
	               'r2'=r2, 
	               'beta'=beta, 
	               'phi'=phi,
	               'n_mods'=n_mods)

	cem@parameters <- c(cem@parameters, params)
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
#' 		  (default) or "signed" networks. Default \code{"unsigned"}.
#' @param cor_function A character string indicating the correlation function 
#' 		  to be used. Default \code{'cor'}.
#' @param cor_method A character string indicating which correlation 
#' 		  coefficient is to be computed. Default \code{"pearson"}
#' @param verbose Logical. If \code{TRUE}, reports analysis steps. Default \code{FALSE}
#'
#' @return A list containing the soft-threshold selected by WGCNA and scale-free model parameters
#' @examples 
#' # Get example expression data
#' data(expr)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr)
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

	expr <- expr_data(cem)
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
#' 
#' @return The phi parameter
#' @examples 
#' # Get example expression data
#' data(expr)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr)
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
#' 		  	  values of R squared to use to calculate the selected beta. 
#'			  Default: 0.1.
#' 
#' @return A vector containing R squared value and the chosen beta parameter.
#' @examples
#' # Get example expression data
#' data(expr)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr)
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
#'
#' @return The value of the network's connectivity.
#' 
#' @examples 
#' # Get example expression data
#' data(expr)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr)
#' # Get modules and beta data
#' cem <- find_modules(cem)
#' # Get network connectivity
#' get_connectivity(cem)
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

#' Calculate co-expression modules
#' 
#' This function takes a \code{CEMiTool} object containing expression values
#' and, together with the given network parameters, returns the given 
#' co-expression modules.
#' @param cem Object of class \code{CEMiTool}.
#' @param beta Selected soft-threshold value
#' @param network_type A character string indicating to use either "unsigned" 
#'        (default) or "signed" networks. Default \code{"unsigned"}.
#' @param tom_type A character string indicating to use either "unsigned" or 
#' 		  "signed" (default) TOM similarity measure.
#' @param cor_function A character string indicating the correlation function 
#'        to be used. Default \code{'cor'}.
#' @param cor_method A character string indicating which correlation 
#'        coefficient is to be computed. Default \code{"pearson"}.
#' 
#' @return Numeric labels assigning genes to modules.
#' 
#' @examples 
#' # Get example expression data
#' data(expr)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr)
#' # Get modules using an example beta value of 7
#' mods <- get_mods(cem, beta=7)
#'
#' @rdname get_mods
#' @export
setGeneric('get_mods', function(cem, ...) {
	standardGeneric('get_mods')
})
#' @rdname get_mods
setMethod('get_mods', signature(cem='CEMiTool'),
	function(cem, beta, network_type="unsigned", tom_type="signed", 
			 cor_function="cor", cor_method="pearson", min_ngen=20) {
    
	expr <- expr_data(cem)
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
	adj <- WGCNA::adjacency(expr_t, power=beta, type=network_type, corFnc=cor_function, corOptions=cor_options)
	cem@adjacency <- adj
			    
	# Calculating Topological Overlap Matrix
	if (tom_type == 'signed') {
		tom <- WGCNA::TOMsimilarity(adj*sign(WGCNA::cor(expr_t)), TOMType=tom_type)
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
#'
#' @return Numeric labels assigning genes to modules.
#'
#' @examples 
#' # Get example expression data
#' data(expr)
#' # Initialize new CEMiTool object with expression data
#' cem <- new_cem(expr)
#' # Get modules using an example beta value of 7
#' mods <- get_mods(cem, beta=7)
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
	    expr <- expr_data(cem)
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
#' @rdname mod_summary
#' @export
setGeneric('mod_summary', function(cem, ...) {
    standardGeneric('mod_summary')
})

#' @rdname mod_summary
setMethod('mod_summary', signature(cem='CEMiTool'),
          function(cem, method=c('mean', 'eigengene'),
                   verbose=FALSE){

              method <- match.arg(method)

		  	  if(length(cem@module) == 0){
				  stop("No modules in CEMiTool object! Did you run find_modules()?")
			  }

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
