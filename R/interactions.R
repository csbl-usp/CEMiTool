#' Includes interaction information to CEMiTool object
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param int_df a data.frame or matrix containing two columns
#'
#' @return Object of class \code{CEMiTool} 
#'
#' @examples
#'
#' @rdname include_interaction
#' @export
setGeneric('include_interaction', function(cem_obj, ...) {
    standardGeneric('include_interaction')
})

#' @rdname find_modules
#' @export
setMethod('include_interaction', signature('CEMiTool'), 
          function(cem_obj, int_df)
{
    genes_by_module <- split(cem_obj@module$genes, cem_obj@module$modules)
    x <- genes_by_module[[1]]
    int <- igraph::graph_from_data_frame(int_df[which(int_df[,1] %in% x | int_df[, 2] %in% x), ])
    lapply(genes_by_module)
    return(cem_obj)
})
