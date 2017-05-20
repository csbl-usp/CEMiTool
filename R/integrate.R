#' Integrates CEMiTool analyses
#'
#' Returns the occurence of edges between different analyses
#'
#' @param analyses List of objects of class \code{CEMiTool}
#'
#' @examples
#' studies <- list('GSE12345'=cem1, 'GSE54321'=cem2)
#' overlaps <- cemoverlap(studies)
#'
#' @export
cemoverlap <- function(analyses) {
    edgelist <- lapply(seq_along(analyses), function(index) {
        cem <- analyses[[index]]
        cem_name <- names(analyses[index])
        # splits by module
        mods <- split(cem@module[,'genes'], cem@module[, 'modules'])
        mods_log <- sapply(mods, length) < 2
        mods <- mods[!mods_log]
        mods <- mods[names(mods) != 'Not.Correlated']
        
        # combines all genes inside each module
        per_mod <- lapply(mods, function(mod) {
               edges <- as.data.table(t(gRbase::combnPrim(mod,2)))
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
    # Sum of cemitool objects containing pair and
    #   order dataframe by sum of occurrences 
    out_presentin <- apply(out[,!colnames(out) %in% c('gene1', 'gene2')], 1, sum)
    out_order <- order(out_presentin, decreasing=TRUE)
    out$presentin <- out_presentin
    out <- out[out_order,]

    return(out)
}
