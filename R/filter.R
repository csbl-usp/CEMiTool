#' @importFrom pracma gammainc
#' @import matrixStats
#'
NULL

#' Filter gene expression table
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param filter_pval P-value cutoff for gene selection.
#' @param n_genes Number of genes to be selected.
#' @param pct Percentage of most expressed genes to keep (before variance filter).
#' @param apply_vst Logical. If TRUE, will apply variance stabilizing transform before filtering data.
#' @param ... Optional parameters.
#'
#' @return Object of class \code{CEMiTool} with selected genes
#'
#' @examples
#' # Get example expression data
#' data(expr)
#' # Initialize CEMiTool object with expression
#' cem <- new_cem(expr)
#' # Filter genes
#' cem1 <- filter_expr(cem)
#' # Check selected genes
#' head(expr_data(cem1))
#' # Filter genes and apply variance stabilizing transformation
#' cem2 <- filter_expr(cem, apply_vst=TRUE)
#' # Check results
#' head(expr_data(cem2))
#'
#' @rdname filter_expr
#' @export
setGeneric('filter_expr', function(cem, ...) {
    standardGeneric('filter_expr')
})

#' @rdname filter_expr
#' @export
setMethod('filter_expr', signature('CEMiTool'),
    function(cem, filter_pval=0.1, n_genes, pct=0.75, apply_vst=FALSE){
        if(!missing(n_genes)){
            if(!missing(filter_pval)){
                stop("Please specify exclusively filter_pval or n_genes.")
            }
        }
        expr <- expr_data(cem, filtered=FALSE)
        if(nrow(expr) == 0){
            stop("CEMiTool object has no expression file!")
        }

        #vars <- mget(ls())
        #vars$expr <- NULL
        #cem <- get_args(cem=cem, vars=vars)

        expr <- expr_pct_filter(expr, pct)

        temp <- as.matrix(expr)
        rownames(temp) <- rownames(expr)
        colnames(temp) <- names(expr)
        expr <- temp

        expr_var <- matrixStats::rowVars(expr)

        expr <- expr[which(expr_var!=0),]

        if(!missing(n_genes)){
            n_genes <- min(n_genes, nrow(expr))
        }

        if (apply_vst){
            expr <- vst(expr)
            temp <- data.frame(expr)
            rownames(temp) <- rownames(expr)
            names(temp) <- colnames(expr)
            expr_data(cem) <- temp
        }

        expr_var <- matrixStats::rowVars(expr)
        names(expr_var) <- rownames(expr)

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
            filter_pval <- sort(p)[n_genes]
            names(filter_pval) <- NULL
        }

        selected <- which(p <= filter_pval)

        if (length(selected) > 0) {
            cem@selected_genes <- names(selected)
        } else {
            cem@selected_genes <- selected
            warning('No gene left after the filtering')
        }

        cem@parameters <- c(cem@parameters,
                            n_genes=length(selected),
                            filter_pval=filter_pval)
        return(cem)
    })

#' Perform variance stabilizing transformation on expression file.
#'
#' @keywords internal
#'
#' @param expr expression file containing genes in the rows and samples in the columns
#'
#' @return  A data.frame containing the results.
vst <- function(expr) {
    #gene_mean <- apply(expr, 1, mean)
    #gene_var  <- apply(expr, 1, var)
    
    gene_mean <- matrixStats::rowMeans2(expr)
    gene_var <- matrixStats::rowVars(expr)

    if(WGCNA::cor(gene_mean, gene_var, method="spearman") > 0.5){
        r <- sum(gene_mean^4)/(sum(gene_var*(gene_mean^2)) - sum(gene_mean^3))
        return(sqrt(r)*asinh(sqrt(expr/r)))
    }else{
        return(expr)
    }
}

#' Filter genes based on expression.
#'
#' @keywords internal
#'
#' @param expr expression file containing genes in the rows and samples in the columns
#' @param pct percentage of most expressed genes to maintain
#'
#' @return A data.frame containing the results
expr_pct_filter <- function(expr, pct=0.75){
    rows <- floor(pct*nrow(expr))
    val <- apply(expr, 1, mean)
    sel_rows <- order(val, decreasing=TRUE)[1:rows]
    expr <- expr[sel_rows, ]
    return(expr)
}
