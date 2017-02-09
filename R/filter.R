#' Filter gene expression table
#'
#' @param cem Object of class \code{CEMiTool}. 
#' @param pval P-value cutoff for gene selection.
#' @param n_genes Number of genes to be selected.
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with selected genes 
#'
#' @example
#  filtered_cem <- filter_expr(cem)
#'
#'
#' @rdname filter_expr
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
    var_var <- var(expr_var)

    ah <- mean_var/var_var + 2*var_var
    bh <- mean_var*(ah - 1)

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
