#' CEMiTool report
#'
#' Creates report for CEMiTool results
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param title Report title.
#' @param directory Directory name for results.
#'
#' @example
#'
#' @rdname generate_report
#' @export
setGeneric('generate_report', function(cem_obj, ...) {
    standardGeneric('generate_report')
})

#' @rdname generate_report
setMethod('generate_report', signature('CEMiTool'),
          function(cem_obj, title="CEMiTool Report", directory="./reports")
{
    # Create the index page
    index_page <- ReportingTools::HTMLReport(shortName = "index",
                             title = title,
                             reportDirectory = directory)
    # Create the module page
    mod_rep <- module_report(cem_obj, directory)
    # add a link to module page in index page
    ReportingTools::publish(ReportingTools::Link(mod_rep, report=index_page), index_page)
    if(!is.null(cem_obj@enrichment)){
        # Create the enrichment page
        enrich_rep <- enrichment_report(cem_obj, directory)
        # add a link to enrichment page in index page
        ReportingTools::publish(ReportingTools::Link(enrich_rep,
                                                     report=index_page), index_page)
    }
    if(nrow(cem_obj@ora)!=0){
        # Create the over representation analysis page
        ora_rep <- ora_report(cem_obj, directory)
        # add a link to ora page in index page
        ReportingTools::publish(ReportingTools::Link(ora_rep,
                                                     report=index_page), index_page)
    }
    ReportingTools::finish(index_page)
    return(index_page)
})


#' Module analysis report
#'
#' Creates report for modules analysis results
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param directory Directory name for results.
#'
#' @example
#'
#'
#' @rdname module_report
#' @export
setGeneric('module_report', function(cem_obj, ...) {
    standardGeneric('module_report')
})


#' @rdname module_report
setMethod('module_report', signature('CEMiTool'),
          function(cem_obj, directory)
{
    mod_rep <- ReportingTools::HTMLReport(shortName = "mod_rep",
                          title="Modules",
                          reportDirectory = directory)
    mod_df <- as.data.frame(table(cem_obj@module$modules))
    colnames(mod_df) <- c("Module", "No.Genes")
    ReportingTools::publish(mod_df, mod_rep)
    ReportingTools::finish(mod_rep)
    return(mod_rep)
})

#' Enrichment analysis report
#'
#' Creates report for modules enrichment analysis results
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param directory Directory name for results.
#'
#' @example
#'
#'
#' @rdname enrichment_report
#' @export
setGeneric('enrichment_report', function(cem_obj, ...) {
    standardGeneric('enrichment_report')
})

#' @rdname enrichment_report
setMethod('enrichment_report', signature('CEMiTool'),
          function(cem_obj, directory)
{
    nes <- cem_obj@enrichment$nes
    pv <- cem_obj@enrichment$pval

    colnames(nes)[-1] <- paste0(colnames(nes)[-1], ": NES")
    colnames(pv)[-1] <- paste0(colnames(pv)[-1], ": p-value")

    # join p-values and nes
    enrich_df <- merge(nes, pv, by.x="pathway")

    # reorder columns
    enrich_df <- enrich_df[, c("pathway", sort(colnames(enrich_df)[-1]))]
    enrich_rep <- ReportingTools::HTMLReport(shortName = "enrich_rep",
                          title="Gene Set Enrichment Analysis",
                          reportDirectory = directory)
    ReportingTools::publish(cem_obj@enrichment_plot, enrich_rep)
    ReportingTools::publish(enrich_df, enrich_rep)
    ReportingTools::finish(enrich_rep)
    return(enrich_rep)
})

#' Over-representation analysis report
#'
#' Creates report for modules over-representation analysis results
#'
#' @param cem_obj Object of class \code{CEMiTool}.
#' @param directory Directory name for results.
#' @param columns Character vector with column names to be shown in report.
#'
#' @example
#'
#' @rdname ora_report
#' @export
setGeneric('ora_report', function(cem_obj, ...) {
    standardGeneric('ora_report')
})

#' @rdname ora_report
setMethod('ora_report', signature('CEMiTool'),
          function(cem_obj, directory,
                   columns=c("ID", "Count", "GeneRatio",
                                 "BgRatio", "p.adjust"))
{
    ora_list <- split(cem_obj@ora[, columns], cem_obj@ora$Module)

    ora_rep <- ReportingTools::HTMLReport(shortName = "ora_rep",
                          title="Over Representation Analysis",
                          reportDirectory = directory)
    for(mod in names(ora_list)){
        ReportingTools::publish(hwriter::hwrite(paste("Module", mod), heading=3), ora_rep)
        if(!is.null(cem_obj$barplot_ora)) {
            ReportingTools::publish(hwriter::hwrite("Barplot", heading=4), ora_rep)
            ReportingTools::publish(cem_obj$barplot_ora[[mod]]$pl, ora_rep)
        }
        ReportingTools::publish(hwriter::hwrite("Table", heading=4), ora_rep)
        row_order <- order(ora_list[[mod]][, "p.adjust"], decreasing = F)
        ReportingTools::publish(ora_list[[mod]][row_order,],
                                ora_rep, name=mod)
    }
    ReportingTools::finish(ora_rep)
    return(ora_rep)
})
