#' Filter gene expression table
#'
#' @param cem Object of class \code{CEMiTool}. 
#' @param pval P-value cutoff for gene selection.
#' @param n_genes Number of genes to be selected.
#' @param ... Other optional parameters.
#'
#' @return Object of class \code{CEMiTool} with selected genes 
#'
#' @rdname filter_expr
#'
#' @example
#' filter_expr(cem) 
#'
#' @export
setGeneric('filter_expr', function(cem, ...) {
    standardGeneric('filter_expr')
})


#' @rdname filter_expr
#' @export
setMethod('filter_expr', signature('CEMiTool'),
          function(cem, pval=0.05, n_genes)
{
    if(!missing(n_genes)){
        if(!missing(pval)){
            stop("Please specify exclusively pval or n_genes.")
        }
    }
    expr <- expr_data(cem, filtered=FALSE)

    if(!missing(n_genes)){
        n_genes <- min(n_genes, nrow(expr))
    }
    expr_var <- apply(expr, 1, var)

    mean_var <- mean(expr_var)
    squared_mean_var <- mean(expr_var)^2
    mean_var_sqr <- mean(expr_var^2)

    ah <- squared_mean_var/(mean_var_sqr - squared_mean_var) + 2
    bh <- (ah-1)*(ah-2)*(mean_var_sqr - squared_mean_var)/mean_var

    p <- sapply(expr_var, function(x) {
        ig <- pracma::gammainc(bh/x, ah)['uppinc']
        g <- gamma(ah)
        return(1 - ig/g)
    })

    names(p) <- gsub('.uppinc', '', names(p))

    if(!missing(n_genes)){
        pval <- sort(p)[n_genes]
        names(pval) <- NULL
    }

    selected <- which(p <= pval)
    
    if (length(selected) > 0) {
        cem@selected_genes <- names(selected)
    } else {
        cem@selected_genes <- selected
        warning('No gene left after the filtering')
    }

    cem@parameters <- c(cem@parameters, 
                            n_genes=length(selected),
                            filter_pval=pval)

    return(cem)
})
