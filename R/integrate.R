#' Integrates CEMiTool analyses
#'
#' Returns the occurence of edges between different analyses
#'
#' @param analyses List of objects of class \code{CEMiTool}
#'
#' @examples
#' studies <- list('GSE12345'=cem_obj1, 'GSE54321'=cem_obj2)
#' overlaps <- cemoverlap(studies)
#'
#' @export
cemoverlap <- function(analyses) {
    edgelist <- lapply(seq_along(analyses), function(index) {
        cem_obj <- analyses[[index]]
        cem_name <- names(analyses[index])
        # splits by module
        mods <- split(cem_obj@module, cem_obj@module[, 'modules'])
        
        # combines all genes inside each module
        per_mod <- lapply(mods, function(mod) {
               edges <- as.data.table(t(gRbase::combnPrim(mod[,1],2)))
               return(edges)
        })
        edges <- do.call(rbind, per_mod)
        setnames(edges, c('gene1', 'gene2'))
        edges[, eval(cem_name) := TRUE]
        return(edges)
    })

    # merges all studies
    out <- Reduce(function(...) { 
                      merge(..., by=c('gene1', 'gene2'),
                                       all=TRUE) },
                 edgelist)
    for(col in colnames(out)) {
        set(out,which(is.na(out[[col]])),col,FALSE)
    }
    setDF(out)
    return(out)
}
