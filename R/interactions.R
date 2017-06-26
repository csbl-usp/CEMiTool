#' @importFrom igraph graph_from_data_frame
NULL

#' Retrieve and set interaction data to CEMiTool object
#'
#' @param cem Object of class \code{CEMiTool}.
#' @param value a data.frame or matrix containing two columns
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
          function(cem) {
              return(cem@interactions)
          })
#' @rdname interactions_data
#' @export
setGeneric("interactions_data<-", function(cem, value) {
    standardGeneric("interactions_data<-")
})

#' @rdname interactions_data
setReplaceMethod("interactions_data", signature("CEMiTool"),
                 function(cem, value){
                     if(nrow(cem@module) == 0){
                         stop("No genes in modules! Did you run find_modules?")
                     }
                     genes_by_module <- split(cem@module$genes, cem@module$modules)
                     cem@interactions <- lapply(genes_by_module, function(x) {
                         rows <- which(value[,1] %in% x | value[, 2] %in% x)
                         ig <- igraph::simplify(igraph::graph_from_data_frame(value[rows,], directed=FALSE))
                         return(ig)
                     })
                     return(cem)
                 })