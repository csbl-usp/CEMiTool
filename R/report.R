generate_report <- function(cemitool_obj, title="CEMiTool Report", directory="./reports") {
    # Create the index page
    index_page <- ReportingTools::HTMLReport(shortName = "index",
                             title = title,
                             reportDirectory = directory)
    # Create the module page
    mod_rep <- module_report(cemitool_obj, directory)
    # add a link to module page in index page
    ReportingTools::publish(ReportingTools::Link(mod_rep, report=index_page), index_page)
    if(!is.null(cemitool_obj$enrichment)){
        # Create the enrichment page
        enrich_rep <- enrichment_report(cemitool_obj, directory)
        # add a link to enrichment page in index page
        ReportingTools::publish(ReportingTools::Link(enrich_rep, report=index_page), index_page)
    }
    if(!is.null(cemitool_obj$ora)){
        # Create the over representation analysis page
        ora_rep <- ora_report(cemitool_obj, directory)
        # add a link to ora page in index page
        ReportingTools::publish(ReportingTools::Link(ora_rep, report=index_page), index_page)
    }
    ReportingTools::finish(index_page)
    return(index_page)
}

module_report <- function(cemitool_obj, directory) {
    mod_rep <- ReportingTools::HTMLReport(shortName = "mod_rep",
                          title="Modules",
                          reportDirectory = directory)
    mod_df <- as.data.frame(table(cemitool_obj$module$modules))
    colnames(mod_df) <- c("Module", "No.Genes")
    ReportingTools::publish(mod_df, mod_rep)
    ReportingTools::finish(mod_rep)
    return(mod_rep)
}

enrichment_report <- function(cemitool_obj, directory) {

    nes <- cemitool_obj$enrichment$nes
    pv <- cemitool_obj$enrichment$pval

    colnames(nes)[-1] <- paste0(colnames(nes)[-1], ": NES")
    colnames(pv)[-1] <- paste0(colnames(pv)[-1], ": p-value")

    # join p-values and nes
    enrich_df <- merge(nes, pv, by.x="pathway")

    # reorder columns
    enrich_df <- enrich_df[, c("pathway", sort(colnames(enrich_df)[-1]))]
    enrich_rep <- ReportingTools::HTMLReport(shortName = "enrich_rep",
                          title="Gene Set Enrichment Analysis",
                          reportDirectory = directory)
    ReportingTools::publish(cemitool_obj$enrichment_plot, enrich_rep)
    ReportingTools::publish(enrich_df, enrich_rep)
    ReportingTools::finish(enrich_rep)
    return(enrich_rep)
}

ora_report <- function(cemitool_obj, directory,
                       columns=c("ID", "Count", "GeneRatio",
                                 "BgRatio", "p.adjust")) {
    ora_list <- split(cemitool_obj$ora[, columns], cemitool_obj$ora$Module)

    ora_rep <- ReportingTools::HTMLReport(shortName = "ora_rep",
                          title="Over Representation Analysis",
                          reportDirectory = directory)
    for(mod in names(ora_list)){
        ReportingTools::publish(hwriter::hwrite(paste("Module", mod), heading=3), ora_rep)
        if(!is.null(cemitool_obj$barplot_ora)) {
            ReportingTools::publish(hwriter::hwrite("Barplot", heading=4), ora_rep)
            ReportingTools::publish(cemitool_obj$barplot_ora[[mod]]$pl, ora_rep)
        }
        ReportingTools::publish(hwriter::hwrite("Table", heading=4), ora_rep)
        row_order <- order(ora_list[[mod]][, "p.adjust"], decreasing = F)
        ReportingTools::publish(ora_list[[mod]][row_order,],
                                ora_rep, name=mod)
    }
    ReportingTools::finish(ora_rep)
    return(ora_rep)
}
