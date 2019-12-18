#' @importFrom pracma gammainc
#' @importFrom matrixStats rowVars rowMeans2
#'
NULL

#' Filter a gene expression data.frame
#'
#' @param expr A data.frame containing expression data
#' @param pct Percentage of most expressed genes to keep.
#' @param apply_vst Logical. If TRUE, will apply variance stabilizing transform
#' before filtering data.
#'
#' @details This function takes a gene expression data.frame and applies a
#' percentage filter (keeps the \code{pct}) most expressed genes). If
#' \code{apply_vst} is TRUE, a variance stabilizing transformation is applied on
#' gene expression values as long as mean and variance values have a Spearman's
#' rho of over 0.5. This transformation is intended to remove this dependence
#' between the parameters. One should then apply the \code{select_genes}
#' function to get significant genes.
#'
#' @return A data.frame containing filtered expression data
#'
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Filter genes
#' expr_f <- filter_genes(expr0)
#' # Check selected genes
#' expr_f[1:5, 1:5]
#' # Filter genes and apply variance stabilizing transformation
#' expr_f2 <- filter_genes(expr0, apply_vst=TRUE)
#' # Check results
#' expr_f2[1:5, 1:5]
#' # Selected genes
#' selected <- select_genes(expr_f2)
#' # Get data.frame with only selected genes
#' expr_s <- expr_f2[selected, ]
#' # Check results
#' expr_s[1:5, 1:5]
#' @rdname filter_expr
#' @export
filter_genes <- function(expr, pct=0.75, apply_vst=FALSE){

    if(nrow(expr) == 0){
        stop("No expression data found.")
    }
    if(pct == 0){
        stop("Argument pct cannot be zero as that would remove all genes.")
    }

    expr <- expr_pct_filter(expr, pct)

    temp <- as.matrix(expr)
    rownames(temp) <- rownames(expr)
    colnames(temp) <- names(expr)
    expr <- temp

    expr_var <- matrixStats::rowVars(expr, na.rm=TRUE)

    expr <- expr[which(expr_var!=0),]

    if (apply_vst){
        expr <- vst(expr)
        # temp <- data.frame(expr)
        # rownames(temp) <- rownames(expr)
        # names(temp) <- colnames(expr)
        # expr_data(cem) <- temp
        # expr <- temp
    }
    return(as.data.frame(expr))
}

#' Select genes based on variance
#'
#' @param expr A data.frame containing expression values
#' @param n_genes (Optional) Number of genes to be selected
#' @param filter_pval P-value cutoff for gene selection
#'
#' @return A vector containing the names of selected genes
#' @export
#'
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Filter genes
#' expr_f <- filter_genes(expr0)
#' # Check selected genes
#' expr_f[1:5, 1:5]
#' # Filter genes and apply variance stabilizing transformation
#' expr_f2 <- filter_genes(expr0, apply_vst=TRUE)
#' # Check results
#' expr_f2[1:5, 1:5]
#' # Selected genes
#' selected <- select_genes(expr_f2)
#' # Get data.frame with only selected genes
#' expr_s <- expr_f2[selected, ]
#' # Check results
#' expr_s[1:5, 1:5]
select_genes <- function(expr, n_genes, filter_pval=0.1){
    expr <- as.matrix(expr)
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
    if(length(selected) == 0) warning("No gene left after filtering!")
    return(names(selected))
}


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
