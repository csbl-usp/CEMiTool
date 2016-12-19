#' Co-expression modules definition
#'
#' Defines co-expression modules
#'
#' @param exprs gene expression \code{"data.frame"}
#' @param cor.method a character string indicating which correlation coefficient is to be computed
#'
#' @return just god knows 
#'
#' @examples
#' coex_mods <- find_modules(exprs=expression.df, cor_method='pearson')
#'
#' @export
find_modules <- function(exprs, cor_method=c('pearson', 'spearman')) {
    #goWGCNA + SplitModules
}


