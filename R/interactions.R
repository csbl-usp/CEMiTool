#' Includes interaction information to CEMiTool object
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param int_df a data.frame or matrix containing two columns
#' @param directed Is the graph directed ? Default is FALSE
#' @param ... parameters for igraph::graph_from_data_frame
#'
#' @return Object of class \code{CEMiTool} 
#'
#' @examples
#'
#' @rdname include_interactions
#' @export
setGeneric('include_interactions', function(cem_obj, ...) {
    standardGeneric('include_interactions')
})

#' @rdname include_interactions
#' @export
setMethod('include_interactions', signature('CEMiTool'), 
          function(cem_obj, int_df, directed=FALSE, ...) {
    genes_by_module <- split(cem_obj@module$genes, cem_obj@module$modules)
    cem_obj@interactions <- lapply(genes_by_module, function(x) {
               rows <- which(int_df[,1] %in% x | int_df[, 2] %in% x)
               ig <- igraph::graph_from_data_frame(int_df[rows,], directed=FALSE)
               return(ig)
    })
    return(cem_obj)
})
