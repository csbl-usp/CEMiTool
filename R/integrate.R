#' @importFrom gRbase combnPrim
#'
NULL

#' Integrates CEMiTool analyses
#'
#' Returns the occurrence of edges between different analyses
#'
#' @param ... Any number of objects of class \code{CEMiTool}
#' @param analyses_list List of objects of class \code{CEMiTool}
#' @param fraction The fraction of objects that an edge pair must be present to be selected (default = 1, accepts values from 0-1)
#'
#' @return Dataframe containing edgelist describing common edges between
#'    the networks defined in module slots from \code{CEMiTool} objects
#'
#' @details Method assumes that all genes inside each module are connected to
#'    every other gene from the same module
#'
#' @examples
#' \dontrun{
#' # Run cemitool twice on expr dataset. In each time, one sample will be removed
#' data(expr)
#' set.seed(10)
#' dset1 <- expr[,-sample(1:ncol(expr), 1)]
#' dset2 <- expr[,-sample(1:ncol(expr), 1)]
#' cem1 <- cemitool(dset1) 
#' cem2 <- cemitool(dset2) 
#' cemoverlap_df <- cemoverlap(cem1, cem2)
#' # Can also run with list: cemoverlap_df <- cemoverlap(list(cem1, cem2))
#'  
#' }
#' @export
cemoverlap <- function(..., analyses_list = NULL, fraction = 1){
    # Several different metrics can be derived from this analysis.
    # See if it is interesting to retrieve the 
    analyses <- c(list(...), analyses_list)
    if(is.null(names(analyses))){
      names(analyses) <- paste0('cem',seq_along(analyses)) 
    }
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

