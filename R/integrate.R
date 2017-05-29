#' @importFrom gRbase combnPrim
#'
NULL

#' Integrates CEMiTool analyses
#'
#' Returns the occurrence of edges between different analyses
#'
#' @param analyses List of objects of class \code{CEMiTool}
#' @param fraction of objects an edge pair must be present \code{CEMiTool}
#'
#'
#' @export
cemoverlap <- function(..., analyses_list = NULL, fraction = 1){
    # Several different metrics can be derived from this analysis.
    # See if it is interesting to retrieve the 
    analyses <- c(list(...), analyses_list)
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
    study_names <- colnames(out)[!colnames(out) %in% c('gene1', 'gene2')] 
    presentin <- apply(out[,study_names], 1, sum)
    out_order <- order(presentin, decreasing=TRUE)
    out$presentin <- presentin
    out <- out[out_order,]
    # Keep only edges present in at least the proportion of cemitool objects specified in 'fraction' variable
    out <- subset(out, presentin >= (fraction * length(study_names)))

    return(out)
}

