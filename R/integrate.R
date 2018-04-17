#' @importFrom gRbase combnPrim
#'
NULL

#' Integrates CEMiTool analyses
#'
#' Returns the occurrence of edges between different analyses
#'
#' @param ... Objects of class \code{CEMiTool}, data.frames or character string containing
#' the path to a file with genes and modules.
#' @param analyses_list List of objects of class \code{CEMiTool}, data.frames or character
#' strings containing the path to files with genes and modules.
#' @param fraction The fraction of objects in which an edge pair must be present to 
#' be selected (default = 1, accepts values from 0-1)
#'
#' @return Object of class \code{data.frame} containing edgelist describing common 
#' edges between the networks defined in module slots from the input objects
#'
#' @details The method assumes that all genes inside each module are connected to
#' every other gene from the same module. 
#'
#' @examples
#' \dontrun{
#' # Run cemitool twice on expr dataset. In each time, one sample will be removed
#' data(expr0)
#' set.seed(10)
#' dset1 <- expr0[,-sample(1:ncol(expr0), 1)]
#' dset2 <- expr0[,-sample(1:ncol(expr0), 1)]
#' cem1 <- cemitool(dset1) 
#' cem2 <- cemitool(dset2) 
#' cemoverlap_df <- cem_overlap(cem1, cem2)
#' # Can also be run with a list: cemoverlap_df <- cemoverlap(list(cem1, cem2))
#' 
#' # Different types of objects can be combined as well, and invalid objects will
#' be removed, with a warning
#' cemoverlap_df <- cem_overlap(cem1, cem2, module_genes(cem), NA, 1, NULL)
#' }
#' @export
cem_overlap <- function(..., analyses_list = NULL, fraction = 1){
    # Several different metrics can be derived from this analysis.
    # See if it is interesting to retrieve the 
    analyses <- c(list(...), analyses_list)
    if(is.null(names(analyses))){
      names(analyses) <- paste0('cem',seq_along(analyses)) 
    }

    analyses <- lapply(analyses, function(x){
        if(is.character(x)){
            if(file.exists(x)){
                data.table::fread(x, data.table=FALSE)
            }
        }else if(class(x) == "CEMiTool"){
            module_genes(x)
        }else if(is.data.frame(x)){
            x
        }else{
            warning("List element ", x, " is not a valid CEMiTool module file, 
                    CEMiTool object or data.frame and will be removed.")
            NULL
        }
    })
    analyses <- Filter(Negate(is.null), analyses)

    edgelist <- lapply(seq_along(analyses), function(index){
        cem <- analyses[[index]]
        cem_name <- names(analyses[index])
        # splits by module
        mods <- split(cem[,'genes'], cem[, 'modules'])
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
    presentin <- rowSums(out[,study_names])
    out_order <- order(presentin, decreasing=TRUE)
    out$presentin <- presentin
    out <- out[out_order,]
    # Keep only edges present in at least the proportion of cemitool objects specified in 'fraction' variable
    out <- subset(out, presentin >= (fraction * length(study_names)))

    return(out)
}

