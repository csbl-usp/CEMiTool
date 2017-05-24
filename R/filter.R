#' @import pracma
#' @import stats
#'
NULL

#' Filter gene expression table
#'
#' @param cem Object of class \code{CEMiTool}. 
#' @param pval P-value cutoff for gene selection.
#' @param n_genes Number of genes to be selected.
#' @param pct .
#' @param dtype Type of experiment. One of \code{'microarray'} or
#'           \code{'rnaseq'}. Default is \code{'microarray'}.
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with selected genes 
#'
#' @examples
#' filtered_cem <- filter_expr(cem)
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
        function(cem, pval=0.1, n_genes, pct=0.75, dtype=c('microarray', 'rnaseq')){
              if(!missing(n_genes)){
                  if(!missing(pval)){
                      stop("Please specify exclusively pval or n_genes.")
                  }
              }
              expr <- expr_data(cem, filtered=FALSE)
              
              expr <- expr_pct_filter(expr, pct)
              
              expr_var <- apply(expr, 1, var)
              
              expr <- expr[which(expr_var!=0),]
              
              if(!missing(n_genes)){
                  n_genes <- min(n_genes, nrow(expr))
              }
              
              dtype <- match.arg(dtype)
              
              if (dtype == 'rnaseq') {
                  message("New filter")
                  expr <- vst(expr)
                  expr_data(cem) <- expr
              }
              
              expr_var <- apply(expr, 1, var)
              
              mean_var <- mean(expr_var)
              var_var <- var(expr_var)
              
              ah <- mean_var^2/var_var + 2
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

vst <- function(expr) {
    gene_mean <- apply(expr, 1, mean)
    gene_var  <- apply(expr, 1, var)
    
    if(stats::cor(gene_mean, gene_var, method="spearman") > 0.5){
        r <- sum(gene_mean^4)/(sum(gene_var*(gene_mean^2)) - sum(gene_mean^3))
        return(sqrt(r)*asinh(sqrt(expr/r)))
    }else{
        return(expr)
    }
}

expr_pct_filter <- function(expr, pct=0.75){
    rows <- floor(pct*nrow(expr))
    val <- apply(expr, 1, mean)
    sel_rows <- order(val, decreasing=TRUE)[1:rows]
    expr <- expr[sel_rows, ]
    return(expr)
}
