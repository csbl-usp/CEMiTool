#' @importFrom igraph graph_from_data_frame
NULL

#' Retrieve and set interaction data to CEMiTool object
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param int_df a data.frame or matrix containing two columns
#' @param directed Is the graph directed ? Default is FALSE
#' @param ... parameters for igraph::graph_from_data_frame
#'
#' @return Object of class \code{CEMiTool} 
#'
#' @examples
#' # Get example CEMiTool object
#' cem <- CEMiTool::cem
#' # Read example interactions data
#' int_df <- read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool"))
#' # Insert interactions data
#' cem <- interactions_data(cem, int_df)
#' # Check interactions data
#' interactions_data(cem)
#' 
#' @rdname interactions_data
#' @export
setGeneric('interactions_data', function(cem, ...) {
    standardGeneric('interactions_data')
})

#' @rdname interactions_data
#' @export
setMethod('interactions_data', signature('CEMiTool'), 
          function(cem, int_df, directed=FALSE, ...) {
              genes_by_module <- split(cem@module$genes, cem@module$modules)
              cem@interactions <- lapply(genes_by_module, function(x) {
                                                 rows <- which(int_df[,1] %in% x | int_df[, 2] %in% x)
                                                 ig <- igraph::simplify(igraph::graph_from_data_frame(int_df[rows,], directed=FALSE))
                                                 return(ig)
              })
              return(cem)
          })
