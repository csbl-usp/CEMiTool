#' @importFrom igraph graph_from_data_frame
NULL

#' Includes interaction information to CEMiTool object
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param int_df a data.frame or matrix containing two columns
#' @param directed Is the graph directed ? Default is FALSE
#' @param ... parameters for igraph::graph_from_data_frame
#'
#' @return Object of class \code{CEMiTool} 
#'
#' @examples
#' \dontrun{
#' int_df <- read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool"))
#' cem_with_interactions <- include_interactions(cem, int_df)
#' }
#' @rdname include_interactions
#' @export
setGeneric('include_interactions', function(cem, ...) {
    standardGeneric('include_interactions')
})

#' @rdname include_interactions
#' @export
setMethod('include_interactions', signature('CEMiTool'), 
          function(cem, int_df, directed=FALSE, ...) {
              genes_by_module <- split(cem@module$genes, cem@module$modules)
              cem@interactions <- lapply(genes_by_module, function(x) {
                                                 rows <- which(int_df[,1] %in% x | int_df[, 2] %in% x)
                                                 ig <- igraph::simplify(igraph::graph_from_data_frame(int_df[rows,], directed=FALSE))
                                                 return(ig)
              })
              return(cem)
          })
