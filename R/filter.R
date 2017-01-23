#' Filter gene expression table
#'
#' @param cem_obj Object of class \code{CEMiTool} 
#' @param pval P-value cutoff for gene selection
#'
#' @return Object of class \code{CEMiTool} with selected genes 
#'
#' @rdname filter_expr
#' @export
setGeneric('filter_expr', function(cem_obj, ...) {
    standardGeneric('filter_expr')
})


#' @rdname filter_expr
#' @export
setMethod('filter_expr', signature('CEMiTool'),
          function(cem_obj, pval=0.05)
{
    expr <- cem_obj@expression
    expr_var <- apply(expr, 1, var)

    var_mean <- mean(expr_var)
    var_mean_sqr <- mean(expr_var)

    ah <- (var_mean^2)/((var_mean_sqr - var_mean^2) + 2)
    bh <- (ah-1)*(ah-2)*(var_mean_sqr - var_mean^2)/var_mean

    p <- sapply(expr_var, function(x) {
        ig <- pracma::gammainc(bh/x, ah)['uppinc']
        g <- gamma(ah)
        return(1 - ig/g)
    })

    names(p) <- gsub('.uppinc', '', names(p))
    selected <- which(p < pval)
    
    if (length(selected) > 0) {
        cem_obj@selected_genes <- names(selected)
    } else {
        cem_obj@selected_genes <- selected
    }
    
    return(cem_obj)
})
